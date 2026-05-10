\\ Sanity tests for compute_m. Run with:  gp -q pari/run_tests.gp
\\ Tests verify the PARI implementation matches values pinned in the Rust
\\ reference (src/depth.rs) and the OEIS draft data for A395913.
\\
\\ Note on syntax: at the GP top level, statements are terminated by newlines
\\ unless wrapped in `{...}` braces. Multi-line calls MUST be braced or they
\\ parse as separate (incomplete) statements per line. Don't remove the braces.

read("pari/gap_depth.gp");

failures = 0;

assert_eq(label, actual, expected) = {
    if(actual == expected,
        printf("PASS: %s\n", label)
    ,
        printf("FAIL: %s\n  expected: %s\n  actual:   %s\n", label, expected, actual);
        failures += 1
    )
};

\\ Test 1: N=100 histogram
{
    h = gap_depth_histogram(100);
    h5 = vector(5, i, if(i <= #h, h[i], 0));
    assert_eq("N=100 histogram", h5, [1, 9, 42, 42, 6])
}

\\ Test 2: N=100 m=0 primes
{
    assert_eq("N=100 m=0 primes", gap_depth_class(100, 0), [2])
}

\\ Test 3: N=100 m=1 primes
{
    assert_eq("N=100 m=1 primes",
        gap_depth_class(100, 1),
        [3, 5, 11, 29, 97, 127, 149, 211, 541])
}

\\ Test 4: N=100 m=4 primes
{
    assert_eq("N=100 m=4 primes",
        gap_depth_class(100, 4),
        [113, 199, 271, 283, 313, 461])
}

\\ Test 5: A395913 first 50 terms (m=3 primes)
{
    a395913_first_50 = [19, 23, 43, 47, 67, 73, 101, 107, 109, 131, 139, 151,
                       163, 173, 179, 181, 193, 227, 229, 233, 241, 263, 269,
                       277, 293, 311, 317, 337, 353, 379, 383, 389, 401, 433,
                       449, 463, 467, 487, 491, 499, 503, 509, 563, 569, 577,
                       593, 599, 601, 613, 619];
    \\ 619 is the 113th prime; N=200 is plenty.
    cls = gap_depth_class(200, 3);
    if(#cls < 50,
        printf("FAIL: A395913 first 50 — only %d m=3 primes found in first 200 primes\n", #cls);
        failures += 1
    ,
        assert_eq("A395913 first 50 terms", cls[1..50], a395913_first_50)
    )
}

\\ Final report
{
    if(failures == 0,
        print("\nAll tests passed.")
    ,
        printf("\n%d test(s) failed.\n", failures);
        quit(1)
    )
}
