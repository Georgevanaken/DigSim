# Rebuild the Pascal reference run used by tests/baselines (needs Free Pascal, e.g. C:\FPC).
# Usage (from DigSim_py):
#   powershell -NoProfile -ExecutionPolicy Bypass -File tools\build_pascal_reference.ps1 -Steps 20
# Everything happens in a temp folder, so nothing lands in the repository or a synced folder.
# 'crt' and 'clrscr' are removed from a COPY of the sources: with redirected stdout they hang the program.
param(
    [int]$Steps = 20,
    [string]$Source = (Join-Path $PSScriptRoot '..\..\DigSim_src\versie 78 effect of salt added'),
    [string]$Fpc = 'C:\FPC\bin\i386-win32\fpc.exe',
    [string]$Work = (Join-Path $env:TEMP 'digsim_ref')
)
$ErrorActionPreference = 'Stop'
if (Test-Path $Work) { Remove-Item $Work -Recurse -Force }
New-Item $Work -ItemType Directory | Out-Null
Copy-Item (Join-Path $Source '*.pas') $Work
Copy-Item (Join-Path $PSScriptRoot '..\tests\data\meal.txt') $Work

foreach ($f in 'digsim_v78.pas', 'insulin.pas') {
    $p = Join-Path $Work $f
    $s = Get-Content $p -Raw
    $s = $s -replace '(?i)\bcrt,\s*', '' -replace '(?im)^\s*clrscr;\s*$', ''
    Set-Content $p $s -NoNewline
}

Push-Location $Work
& $Fpc -Mobjfpc -O2 -oreference digsim_v78.pas | Select-String 'Error|Fatal'
if (Test-Path 'reference') { Move-Item 'reference' 'reference.exe' -Force }
Set-Content "$Work\in.txt" "$Steps"
$run = Join-Path $Work "n$Steps"
New-Item $run -ItemType Directory | Out-Null
Copy-Item meal.txt $run
$p = Start-Process "$Work\reference.exe" -WorkingDirectory $run -PassThru -NoNewWindow `
    -RedirectStandardInput "$Work\in.txt" -RedirectStandardOutput "$run\stdout.log"
if (-not $p.WaitForExit(900000)) { $p.Kill(); throw 'Pascal run timed out' }
Pop-Location
Write-Host "Results: $run\Results.txt ($((Get-Content "$run\Results.txt").Count) lines)"
