#!/usr/bin/env pwsh
# Run all statistics for 10^1, ... 10^n

param(
    [Parameter(Mandatory = $false)]
    [ValidateRange(1, [int]::MaxValue)]
    [int]$pow,
    
    [Parameter(Mandatory = $false)]
    [Alias("o")]
    [string]$OutputDir = "out",
    
    [Parameter(Mandatory = $false)]
    [Alias("s")]
    [Alias("--seed-file")]
    [Alias("--seed-set")]
    [string]$SeedSet
)

function runstatistics {
    param (
        $count,
        $DIR
    )
    Write-Host "Running statistics for count = $count"

    # Build base command arguments
    $baseArgs = @("-o", "$DIR")
    $baseArgs += @("--preserve-order")
    if ($SeedSet) {
        $baseArgs += @("--seed-file", $SeedSet)
    } else {
        $baseArgs += @("--generator", "primes")
    }
    if ($count) {
        $baseArgs += @("-n", $count)
    }

    ./target/release/pgd                 @baseArgs
    # ./target/release/pgd locality        @baseArgs
    ./target/release/pgd oeis-export     @baseArgs
    ./target/release/pgd mod-residue 30  @baseArgs
    ./target/release/pgd mod-residue 210 @baseArgs
    ./target/release/pgd class-quantiles @baseArgs
    ./target/release/pgd overlay         @baseArgs
    ./target/release/pgd gap-address     @baseArgs
}

if ($pow) {
    foreach ($x in 1..$pow)
    {
        Write-Host "Running statistics for 10^$x"
        $count = [math]::Pow(10, $x);

        $DIR="$OutputDir/10^$x"
        runstatistics -count $count -DIR $DIR
    }
} else {
    # run once on entire seed file, save directly to out/dir
    runstatistics -count $null -DIR $OutputDir
}