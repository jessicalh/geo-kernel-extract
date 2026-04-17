# launch_reader.ps1 — build and run h5reader on Windows.
#
# Usage (from x64 Native Tools prompt or a PowerShell with VS env set):
#   .\launch_reader.ps1                                     # build + run, no H5
#   .\launch_reader.ps1 path\to\protein_analysis.h5         # with an H5
#
# Override preset via $env:H5READER_PRESET (default: win-msvc).

$ErrorActionPreference = 'Stop'

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -Path $scriptDir

$preset = $env:H5READER_PRESET
if (-not $preset) { $preset = 'win-msvc' }

$buildDir = Join-Path 'build' $preset
if (-not (Test-Path $buildDir)) {
    cmake --preset $preset
}

cmake --build --preset $preset

$bin = Join-Path $buildDir 'h5reader.exe'
if (-not (Test-Path $bin)) {
    Write-Error "h5reader binary missing after build: $bin"
    exit 1
}

& $bin $args
exit $LASTEXITCODE
