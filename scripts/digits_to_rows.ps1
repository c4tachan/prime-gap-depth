#!/usr/bin/env pwsh
# Convert decimal digits from a file to one digit per row
# Usage: ./digits_to_rows.ps1 -InputFile pi_digits.txt -OutputFile seed.txt

param(
    [Parameter(Mandatory = $true)]
    [string]$InputFile,
    
    [Parameter(Mandatory = $false)]
    [string]$OutputFile = "seed_digits.txt"
)

# Read the file content
$content = Get-Content -Path $InputFile -Raw

# Remove decimal point and all non-digit characters (spaces, newlines, etc.)
$digits = $content -replace '[^0-9]', ''

# Split each digit into separate lines
$digitArray = $digits -split '' | Where-Object { $_ -ne '' }

# Write to output file
$digitArray | Set-Content -Path $OutputFile

Write-Host "Extracted $($digitArray.Count) digits from '$InputFile' to '$OutputFile'"
