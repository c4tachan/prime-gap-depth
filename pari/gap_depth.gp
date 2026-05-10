\\ Prime Gap Depth — PARI/GP reference implementation.
\\
\\ See ../docs/algorithm.md for the formal construction. In summary: starting
\\ with a strictly increasing sequence v as a single row, iteratively assign
\\ each row's smallest element (its leader) the current iteration index, then
\\ partition the remaining elements by within-row consecutive gap into new
\\ rows for the next iteration. m(v[i]) is the iteration at which v[i]
\\ becomes a leader.
\\
\\ This file is intended to be `read()` from other GP scripts. For a
\\ self-contained version suitable for OEIS PROG entries see
\\ pari/oeis_prog_snippet.gp.

\\ compute_m(v): given a strictly increasing vector v, return a vector m of
\\ the same length where m[i] is the gap-depth of v[i].
compute_m(v) = {
    my(n = #v);
    if(n == 0, return([]));
    my(m = vector(n, i, -1));
    my(stk = List([[0, vector(n, i, i)]]));
    while(#stk > 0,
        my(top = stk[#stk]); listpop(stk);
        my(lvl = top[1], row = top[2]);
        if(#row == 0, next);
        m[row[1]] = lvl;
        if(#row == 1, next);
        \\ Bucket non-leaders by within-row consecutive gap.
        \\ `seen` preserves first-encountered-gap order so traversal is
        \\ deterministic; `bg` maps gap value -> List of positions.
        my(seen = List(), bg = Map(), cur);
        for(j = 2, #row,
            my(g = v[row[j]] - v[row[j-1]]);
            if(mapisdefined(bg, g, &cur),
                listput(cur, row[j]);
                mapput(bg, g, cur)
            ,
                listput(seen, g);
                mapput(bg, g, List([row[j]]))
            )
        );
        for(k = 1, #seen,
            listput(stk, [lvl + 1, Vec(mapget(bg, seen[k]))])
        )
    );
    m
};

\\ gap_depth_class(N, target): return the primes among the first N primes
\\ whose gap-depth equals `target`, in ascending order.
gap_depth_class(N, target) = {
    my(p = primes(N), m = compute_m(p), out = List());
    for(i = 1, N, if(m[i] == target, listput(out, p[i])));
    Vec(out)
};

\\ gap_depth_histogram(N): return a vector h such that h[k+1] = #{i : m(p_i) = k}
\\ for the first N primes (h is sized to max-depth+1).
gap_depth_histogram(N) = {
    my(p = primes(N), m = compute_m(p), maxm = vecmax(m));
    my(h = vector(maxm + 1, i, 0));
    for(i = 1, N, h[m[i] + 1] += 1);
    h
};
