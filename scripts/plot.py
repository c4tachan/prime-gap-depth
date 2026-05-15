#!/usr/bin/env python3
"""
Interactive CDF plotter for prime gap-depth m-classes.

Usage:
    # TSV from pgd plot-prep:
    python3 scripts/plot.py plot_data.tsv

    # One or more OEIS text files (oeis_m*.txt from pgd oeis-export):
    python3 scripts/plot.py oeis_m0.txt oeis_m1.txt oeis_m2.txt ...

The OEIS files have k and p(k) only; the n-axis and pi(k) y-axis options
are not available in that mode.
"""

import re
import sys
import argparse
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_oeis_files(paths):
    """Load one or more oeis_m*.txt files.  Returns same dict as load_data.
    'pi' is set to an empty array because the OEIS format has no prime_index."""
    data = {}
    for path in paths:
        m_match = re.search(r'oeis_m(\d+)', Path(path).name)
        if m_match is None:
            print(f'Warning: cannot determine m from filename {path!r}, skipping')
            continue
        m = int(m_match.group(1))
        ks, ps = [], []
        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                ks.append(int(parts[0]))
                ps.append(int(parts[1]))
        data[m] = {
            'k':    np.array(ks, dtype=np.float64),
            'pi':   np.array([], dtype=np.float64),  # not available in OEIS format
            'p':    np.array(ps, dtype=np.float64),
            'size': len(ks),
        }
    return data


def load_data(path):
    """Return dict: m -> {k, pi, p, size}  (all numpy float64 arrays)."""
    data = {}
    with open(path, encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        col = {h: i for i, h in enumerate(header)}
        for line in f:
            row = line.strip().split('\t')
            if len(row) < 5:
                continue
            m          = int(row[col['m']])
            approx_k   = int(row[col['rank']])
            class_size = int(row[col['class_size']])
            pi_k       = int(row[col['prime_index']])
            p_k        = int(row[col['prime_value']])
            if m not in data:
                data[m] = {'k': [], 'pi': [], 'p': [], 'size': class_size}
            data[m]['k'].append(approx_k)
            data[m]['pi'].append(pi_k)
            data[m]['p'].append(p_k)

    for m in data:
        data[m]['k']  = np.array(data[m]['k'],  dtype=np.float64)
        data[m]['pi'] = np.array(data[m]['pi'], dtype=np.float64)
        data[m]['p']  = np.array(data[m]['p'],  dtype=np.float64)
    return data

# ---------------------------------------------------------------------------
# Curve fitting
# ---------------------------------------------------------------------------

FITS = ['None', 'Power law  k^a', 'PNT  C·k^a·ln k', 'Linear  a·k+b', 'Quadratic ln']

def _stats(y_actual, y_pred, log_space):
    """
    Compute fit statistics.
    If log_space=True, residuals are in ln(y) — RMSE is in nats, typical_ratio = exp(RMSE).
    Returns dict with keys: n, r2, rmse, typical_ratio (None if linear), max_abs_resid.
    """
    if log_space:
        res   = np.log(y_actual) - np.log(y_pred)
        ss_res = np.sum(res ** 2)
        ss_tot = np.sum((np.log(y_actual) - np.log(y_actual).mean()) ** 2)
    else:
        res    = y_actual - y_pred
        ss_res = np.sum(res ** 2)
        ss_tot = np.sum((y_actual - y_actual.mean()) ** 2)

    n    = len(y_actual)
    r2   = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    rmse = np.sqrt(ss_res / n)
    return {
        'n':             n,
        'r2':            r2,
        'rmse':          rmse,
        'typical_ratio': float(np.exp(rmse)) if log_space else None,
        'max_abs_resid': float(np.max(np.abs(res))),
        'log_space':     log_space,
    }

def fit_curve(x_raw, y_raw, fit_type, n_pts=400):
    """
    Fit a curve to (x_raw, y_raw).
    Returns (x_dense, y_dense, label, stats_dict) or None.
    Fitting is done on log-transformed data where appropriate.
    Small-k region (k < 10) is excluded to avoid transient noise.
    """
    if fit_type == 'None':
        return None

    mask = (x_raw > 10) & (y_raw > 1) & np.isfinite(x_raw) & np.isfinite(y_raw)
    if mask.sum() < 4:
        return None
    x, y = x_raw[mask], y_raw[mask]

    x_lo = max(x.min(), 1.0)
    x_hi = x.max()
    x_dense = np.logspace(np.log10(x_lo), np.log10(x_hi), n_pts)

    try:
        if fit_type == 'Power law  k^a':
            a, c    = np.polyfit(np.log(x), np.log(y), 1)
            y_dense = np.exp(c) * x_dense ** a
            y_pred  = np.exp(c) * x ** a
            label   = f'C·k^{a:.4f}'
            st      = _stats(y, y_pred, log_space=True)

        elif fit_type == 'PNT  C·k^a·ln k':
            lnx   = np.log(x)
            lnlnx = np.log(np.maximum(lnx, 1e-10))
            lny   = np.log(y)
            X     = np.column_stack([lnx, np.ones_like(lnx)])
            coeffs, *_ = np.linalg.lstsq(X, lny - lnlnx, rcond=None)
            a, c  = coeffs
            y_dense = np.exp(c) * (x_dense ** a) * np.log(np.maximum(x_dense, 1.0))
            y_pred  = np.exp(c) * (x ** a) * np.log(np.maximum(x, 1.0))
            label   = f'C·k^{a:.4f}·ln k'
            st      = _stats(y, y_pred, log_space=True)

        elif fit_type == 'Linear  a·k+b':
            a, b    = np.polyfit(x, y, 1)
            y_dense = a * x_dense + b
            y_pred  = a * x + b
            label   = f'{a:.3g}·k + {b:.3g}'
            st      = _stats(y, y_pred, log_space=False)

        elif fit_type == 'Quadratic ln':
            lnx2    = np.log(x) ** 2
            a, b    = np.polyfit(lnx2, y, 1)
            y_dense = a * np.log(x_dense) ** 2 + b
            y_pred  = a * np.log(x) ** 2 + b
            label   = f'{a:.3g}·(ln k)² + {b:.3g}'
            st      = _stats(y, y_pred, log_space=False)

        else:
            return None

        return x_dense, y_dense, label, st

    except Exception:
        return None

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------

COLORS = {
    0: '#888888',   # m=0  grey
    1: '#9467bd',   # m=1  purple
    2: '#1f77b4',   # m=2  blue
    3: '#2ca02c',   # m=3  green
    4: '#ff7f0e',   # m=4  orange
    5: '#d62728',   # m=5  red
    6: '#8c564b',   # m=6  brown
    7: '#e377c2',   # m=7  pink
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Interactive m-class CDF plotter for prime gap-depth data')
    parser.add_argument('inputs', nargs='+',
                        help='plot_data.tsv  OR  one or more oeis_m*.txt files')
    args = parser.parse_args()

    # Detect format: single .tsv → TSV mode; anything else → OEIS mode
    if len(args.inputs) == 1 and args.inputs[0].lower().endswith('.tsv'):
        print(f'Loading {args.inputs[0]} ...')
        data = load_data(args.inputs[0])
        oeis_mode = False
    else:
        print(f'Loading {len(args.inputs)} OEIS file(s) ...')
        data = load_oeis_files(args.inputs)
        oeis_mode = True

    m_values = sorted(data.keys())
    print(f'  Classes: {m_values}')
    for m in m_values:
        d = data[m]
        pi_note = ''  if not oeis_mode else '  (no prime_index)'
        print(f'  m={m}: {len(d["k"])} samples, class_size={d["size"]:,}{pi_note}')

    # -------------------------------------------------------------------------
    # Figure layout
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(15, 10))
    fig.patch.set_facecolor('#f8f8f8')

    # Main plot area (right 70%, upper ~67%)
    ax = fig.add_axes([0.30, 0.33, 0.67, 0.62])
    ax.set_facecolor('#ffffff')

    # Stats table area (right 70%, lower ~25%)
    ax_stats = fig.add_axes([0.30, 0.03, 0.67, 0.26])
    ax_stats.set_facecolor('#f0f0f0')
    ax_stats.axis('off')

    # ── control widgets (left 28%) ────────────────────────────────────────────
    ctrl_l, ctrl_w = 0.01, 0.27

    # m-class checkboxes
    check_labels  = [f'm={m}  ({data[m]["size"]:,})' for m in m_values]
    check_default = [m >= 2 for m in m_values]
    ax_check = fig.add_axes([ctrl_l, 0.60, ctrl_w, 0.35])
    ax_check.set_title('m-classes', fontsize=9, pad=2)
    check = mwidgets.CheckButtons(ax_check, check_labels, check_default)
    # Smaller tick font
    for lbl in check.labels:
        lbl.set_fontsize(8)

    # X-axis selection
    x_opts = ['k  — rank within m-class', 'n  — prime index (all primes)']
    ax_xsel = fig.add_axes([ctrl_l, 0.50, ctrl_w, 0.09])
    ax_xsel.set_title('X-axis', fontsize=9, pad=2)
    radio_x = mwidgets.RadioButtons(ax_xsel, x_opts, active=0)
    for lbl in radio_x.labels:
        lbl.set_fontsize(8)

    # Y-axis selection
    y_opts = ['prime value  p(k)', 'prime index  \u03c0(k)']
    ax_ysel = fig.add_axes([ctrl_l, 0.41, ctrl_w, 0.09])
    ax_ysel.set_title('Y-axis', fontsize=9, pad=2)
    radio_y = mwidgets.RadioButtons(ax_ysel, y_opts, active=0)
    for lbl in radio_y.labels:
        lbl.set_fontsize(8)

    # Axis scale
    scale_opts = ['log-log', 'log x / lin y', 'lin x / log y', 'lin-lin']
    ax_scale = fig.add_axes([ctrl_l, 0.22, ctrl_w, 0.18])
    ax_scale.set_title('Axis scale', fontsize=9, pad=2)
    radio_scale = mwidgets.RadioButtons(ax_scale, scale_opts, active=0)
    for lbl in radio_scale.labels:
        lbl.set_fontsize(8)

    # Curve fit
    ax_fit = fig.add_axes([ctrl_l, 0.02, ctrl_w, 0.19])
    ax_fit.set_title('Curve fit', fontsize=9, pad=2)
    radio_fit = mwidgets.RadioButtons(ax_fit, FITS, active=0)
    for lbl in radio_fit.labels:
        lbl.set_fontsize(8)

    # -------------------------------------------------------------------------
    # Update / redraw
    # -------------------------------------------------------------------------
    def update(_=None):
        ax.cla()
        ax_stats.cla()
        ax_stats.axis('off')

        x_choice     = radio_x.value_selected
        y_choice     = radio_y.value_selected
        scale_choice = radio_scale.value_selected
        fit_choice   = radio_fit.value_selected

        use_n_axis   = x_choice.startswith('n')
        use_pi_y     = 'prime index' in y_choice
        unavailable  = (use_n_axis or use_pi_y) and oeis_mode

        # Axis scales
        xlog = scale_choice.startswith('log')
        ylog = 'log y' in scale_choice or scale_choice == 'log-log'
        ax.set_xscale('log' if xlog else 'linear')
        ax.set_yscale('log' if ylog else 'linear')

        xlabel = 'n  — prime index (rank among all primes)' if use_n_axis else 'k  — rank within m-class'
        ylabel = 'p(k)  (prime value)' if 'prime value' in y_choice else '\u03c0(k)  (prime index)'
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)

        if unavailable:
            ax.text(0.5, 0.5,
                    'Not available in OEIS mode\n(no prime_index column)',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=13, color='#888888')
            fig.canvas.draw_idle()
            return
        ax.grid(True, which='both', alpha=0.25, linewidth=0.6)
        ax.grid(True, which='major', alpha=0.5, linewidth=0.8)

        stats_rows = []   # list of dicts, one per active m-class with a fit

        any_plotted = False
        for idx, m in enumerate(m_values):
            if not check.get_status()[idx]:
                continue
            d     = data[m]
            x     = d['pi'] if use_n_axis else d['k']
            y_arr = d['p'] if 'prime value' in y_choice else d['pi']
            if len(x) == 0 or len(y_arr) == 0:
                continue
            col   = COLORS.get(m, f'C{m}')
            n     = data[m]['size']

            ax.scatter(x, y_arr, s=5, alpha=0.55, color=col, zorder=3,
                       label=f'm={m}  (n={n:,})')
            any_plotted = True

            # Fitted curve overlay
            result = fit_curve(x, y_arr, fit_choice)
            if result is not None:
                xf, yf, flabel, st = result
                ax.plot(xf, yf, '-', color=col, linewidth=1.8, alpha=0.9,
                        zorder=4, label=f'  fit: {flabel}')
                stats_rows.append({'m': m, 'color': col, 'label': flabel, 'st': st})

        if any_plotted:
            ax.legend(fontsize=8, loc='upper left', framealpha=0.9,
                      markerscale=2.5, handlelength=1.5)
        else:
            ax.text(0.5, 0.5, 'No classes selected', transform=ax.transAxes,
                    ha='center', va='center', fontsize=13, color='#888888')

        # ── stats table ───────────────────────────────────────────────────────
        if stats_rows:
            log_space = stats_rows[0]['st']['log_space']
            rmse_label = 'RMSE (ln)' if log_space else 'RMSE'

            # Header
            headers = ['m', 'formula', 'n fit', 'R²', rmse_label,
                       'typ. ratio' if log_space else 'max |err|',
                       'max |resid|']
            col_xs  = [0.00, 0.06, 0.28, 0.40, 0.52, 0.65, 0.82]
            col_aligns = ['left', 'left', 'right', 'right', 'right', 'right', 'right']

            y_top   = 0.93
            row_h   = 0.13
            hdr_fs  = 8
            row_fs  = 8.5

            # Draw header row
            for hdr, cx, align in zip(headers, col_xs, col_aligns):
                ax_stats.text(cx, y_top, hdr, transform=ax_stats.transAxes,
                              fontsize=hdr_fs, fontweight='bold',
                              ha=align, va='top', family='monospace',
                              color='#333333')

            # Separator line
            ax_stats.plot([0.0, 1.0], [y_top - 0.02, y_top - 0.02],
                          color='#aaaaaa', linewidth=0.8,
                          transform=ax_stats.transAxes, clip_on=False)

            for i, row in enumerate(stats_rows):
                st  = row['st']
                col = row['color']
                y_r = y_top - row_h * (i + 1)

                n_fit  = f"{st['n']:,}"
                r2_s   = f"{st['r2']:.5f}" if np.isfinite(st['r2']) else '—'
                rmse_s = f"{st['rmse']:.4f}" if log_space else f"{st['rmse']:.3g}"
                if log_space:
                    extra_s = f"{st['typical_ratio']:.4f}×"
                else:
                    extra_s = f"{st['max_abs_resid']:.3g}"
                maxr_s = f"{st['max_abs_resid']:.4f}" if log_space else f"{st['max_abs_resid']:.3g}"

                cells = [f"m={row['m']}", row['label'], n_fit, r2_s,
                         rmse_s, extra_s, maxr_s]

                for cell, cx, align in zip(cells, col_xs, col_aligns):
                    ax_stats.text(cx, y_r, cell, transform=ax_stats.transAxes,
                                  fontsize=row_fs, ha=align, va='top',
                                  family='monospace', color=col)

            # Footnote
            if log_space:
                note = ('R² and RMSE computed in ln(y) space.  '
                        'Typical ratio = exp(RMSE) — geometric mean absolute factor error.')
            else:
                note = 'R² and RMSE computed in y space.'
            ax_stats.text(0.0, 0.01, note, transform=ax_stats.transAxes,
                          fontsize=7, color='#666666', va='bottom')
        else:
            if fit_choice != 'None' and any_plotted:
                ax_stats.text(0.5, 0.5, 'Insufficient data for fit statistics.',
                              transform=ax_stats.transAxes,
                              ha='center', va='center', fontsize=10, color='#888888')
            elif fit_choice == 'None' and any_plotted:
                ax_stats.text(0.5, 0.5, 'Select a curve fit to see statistics.',
                              transform=ax_stats.transAxes,
                              ha='center', va='center', fontsize=10, color='#aaaaaa')

        fig.canvas.draw_idle()

    # Wire up all controls
    check.on_clicked(update)
    radio_x.on_clicked(update)
    radio_y.on_clicked(update)
    radio_scale.on_clicked(update)
    radio_fit.on_clicked(update)

    update()  # initial draw
    plt.show()


if __name__ == '__main__':
    main()
