\\ Emit `pi,p,m` for the first N primes, one per line, comma-separated.
\\ Used by tests/cross_check.sh to diff against the Rust binary's results.csv.
\\ N defaults to 1000. To override:
\\   echo 'N=5000; read("pari/cross_check.gp")' | gp -q

read("pari/gap_depth.gp");

\\ If N is undefined, GP treats it as a polynomial variable (t_POL); use that
\\ as the "not set" sentinel and fall back to the default.
if(type(N) == "t_POL", N = 1000);

{
    p = primes(N);
    m = compute_m(p);
    for(i = 1, N, printf("%d,%d,%d\n", i, p[i], m[i]))
}
