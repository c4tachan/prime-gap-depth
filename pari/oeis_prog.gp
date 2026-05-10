\\ Self-contained PARI/GP program for OEIS gap-depth sequences.
\\ This file IS the PROG snippet you paste into the OEIS submission, minus
\\ the `(PARI)` language tag and `\\ ~~~~` signature line which OEIS expects.
\\
\\ The function `gap_depth_class(N, target)` returns the primes among the
\\ first N primes whose gap-depth equals `target`. Adjust the final line
\\ for each sequence.

gap_depth_class(N, target) = {
    my(p = primes(N), n = #p, m = vector(n, i, -1));
    my(stk = List([[0, vector(n, i, i)]]));
    while(#stk > 0,
        my(top = stk[#stk]); listpop(stk);
        my(lvl = top[1], row = top[2]);
        if(#row == 0, next);
        m[row[1]] = lvl;
        if(#row == 1, next);
        my(seen = List(), bg = Map(), cur);
        for(j = 2, #row,
            my(g = p[row[j]] - p[row[j-1]]);
            if(mapisdefined(bg, g, &cur),
                listput(cur, row[j]); mapput(bg, g, cur)
            ,
                listput(seen, g); mapput(bg, g, List([row[j]]))
            )
        );
        for(k = 1, #seen,
            listput(stk, [lvl + 1, Vec(mapget(bg, seen[k]))])
        )
    );
    my(out = List());
    for(i = 1, n, if(m[i] == target, listput(out, p[i])));
    Vec(out)
};

\\ A395913: primes of gap-depth 3, first 50 terms (619 is the 113th prime).
print(gap_depth_class(200, 3)[1..50])
