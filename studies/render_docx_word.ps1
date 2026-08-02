param(
    [Parameter(Mandatory = $true)]
    [string]$InputDocx,
    [Parameter(Mandatory = $true)]
    [string]$OutputPdf,
    [Parameter(Mandatory = $true)]
    [string]$StatusFile
)

$ErrorActionPreference = "Stop"
$word = $null
$document = $null
try {
    $word = New-Object -ComObject Word.Application
    $word.Visible = $false
    $word.DisplayAlerts = 0
    $word.AutomationSecurity = 3
    $document = $word.Documents.Open($InputDocx, $false, $true, $false)
    $document.SaveAs2($OutputPdf, 17)
    "OK`n$OutputPdf" | Set-Content -LiteralPath $StatusFile -Encoding UTF8
}
catch {
    "ERROR`n$($_.Exception.ToString())" | Set-Content -LiteralPath $StatusFile -Encoding UTF8
    exit 1
}
finally {
    if ($document -ne $null) {
        $document.Close($false)
    }
    if ($word -ne $null) {
        $word.Quit()
    }
}
