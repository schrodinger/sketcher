# Runs a batch file (e.g. vcvarsall.bat or setvars.bat) in cmd.exe and writes
# any environment variables that it adds or changes to $GITHUB_ENV so that they
# are visible to subsequent steps of the job
param(
    [Parameter(Mandatory)] [string] $BatchFile,
    [string[]] $Arguments = @()
)

$ErrorActionPreference = 'Stop'

if (-not (Test-Path $BatchFile)) { throw "$BatchFile not found" }
$batchName = Split-Path -Leaf $BatchFile

$envBefore = New-Object 'System.Collections.Generic.Dictionary[String,String]' ([StringComparer]::OrdinalIgnoreCase)
Get-ChildItem env: | ForEach-Object { $envBefore[$_.Name] = $_.Value }

$output = cmd /c "`"$BatchFile`" $($Arguments -join ' ') > nul && set"
if ($LASTEXITCODE -ne 0) { throw "$batchName failed with exit code $LASTEXITCODE" }

foreach ($line in $output) {
    if ($line -notmatch '^([^=]+)=(.*)$') { continue }
    $name = $matches[1]
    $value = $matches[2]

    # Skip variables that contain a multi-line value (e.g. VERSIONS_JSON)
    # since the batch file shouldn't modify them
    if ($envBefore.ContainsKey($name) -and $envBefore[$name].Contains("`n")) {
        continue
    }

    if (-not $envBefore.ContainsKey($name) -or $envBefore[$name] -ne $value) {
        Add-Content -Path $env:GITHUB_ENV -Value "$name=$value"
    }
}
