# Build and run the Windows Reader. Pass -Package as the first argument to
# build the NSIS installer instead.

$ErrorActionPreference = 'Stop'

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $scriptDir

function RequiredPath([string]$environmentName, [string]$defaultPath) {
    $path = [Environment]::GetEnvironmentVariable($environmentName)
    if (-not $path) { $path = $defaultPath }
    if (-not (Test-Path -LiteralPath $path)) {
        throw "$environmentName is required and does not exist: $path"
    }
    return (Resolve-Path -LiteralPath $path).Path.Replace('\', '/')
}

$cmake = RequiredPath 'H5READER_CMAKE' 'C:\Qt\Tools\CMake_64\bin\cmake.exe'
$qtRoot = RequiredPath 'H5READER_QT_DIR' 'C:\Qt\6.10.2\msvc2022_64'
$vtkRoot = RequiredPath 'H5READER_VTK_DIR' 'C:\Projects\VTK'
$vcpkgRoot = RequiredPath 'VCPKG_ROOT' 'C:\Projects\vcpkg'
$modelRoot = RequiredPath `
    'H5READER_EXPERIMENTAL_SHIELDING_ML_ROOT' `
    'C:\Projects\reader-data\experimental-shielding-ml-bundle-20260722-F006-R007-v1'

if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue)) {
    $vswhere = RequiredPath `
        'H5READER_VSWHERE' `
        'C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe'
    $visualStudio = & $vswhere -latest -products * `
        -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
        -property installationPath
    if (-not $visualStudio) {
        throw 'A Visual Studio installation with the C++ toolchain is required.'
    }
    $devShell = Join-Path $visualStudio 'Common7\Tools\Microsoft.VisualStudio.DevShell.dll'
    Import-Module $devShell
    Enter-VsDevShell -VsInstallPath $visualStudio `
        -DevCmdArguments '-arch=x64' -SkipAutomaticLocation
}

$preset = $env:H5READER_PRESET
if (-not $preset) { $preset = 'win-rwdi' }
$buildDir = Join-Path 'build' $preset
$toolchain = "$vcpkgRoot/scripts/buildsystems/vcpkg.cmake"

& $cmake --preset $preset `
    "-DCMAKE_TOOLCHAIN_FILE=$toolchain" `
    "-DH5READER_QT_DIR=$qtRoot" `
    "-DH5READER_VTK_DIR=$vtkRoot" `
    "-DH5READER_EXPERIMENTAL_SHIELDING_ML_ROOT=$modelRoot"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

& $cmake --build --preset $preset
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

$readerArguments = @($args)
if ($readerArguments.Count -gt 0 -and $readerArguments[0] -eq '-Package') {
    $cpack = Join-Path (Split-Path -Parent $cmake) 'cpack.exe'
    & $cpack --config (Join-Path $buildDir 'CPackConfig.cmake') -G NSIS
    exit $LASTEXITCODE
}

$env:PATH = @(
    (Join-Path $qtRoot 'bin')
    (Join-Path $vtkRoot 'bin')
    (Join-Path $vcpkgRoot 'installed\x64-windows\bin')
    $env:PATH
) -join [IO.Path]::PathSeparator

$binary = Join-Path $buildDir 'h5reader.exe'
if (-not (Test-Path -LiteralPath $binary)) {
    throw "h5reader binary is missing after build: $binary"
}

& $binary $readerArguments
exit $LASTEXITCODE
