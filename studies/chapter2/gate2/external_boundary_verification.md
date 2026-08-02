# Gate 2E: external-boundary verification

**Status:** PASS
**Checks:** 12/12
**Policy:** `observed_interchange`

## Interface evidence

| interface   | policy               |   checks_passed |   checks_total |   maximum_source_series_error_mw_or_pu | observed_hours_outside_ttc   |   source_label_anomaly_events |   maximum_observed_import_mw |   maximum_observed_export_mw | maximum_import_ttc_mw   | maximum_export_ttc_mw   |
|:------------|:---------------------|----------------:|---------------:|---------------------------------------:|:-----------------------------|------------------------------:|-----------------------------:|-----------------------------:|:------------------------|:------------------------|
| AB          | observed_interchange |               6 |              6 |                                      0 |                              |                             2 |                          555 |                          744 |                         |                         |
| US          | observed_interchange |               6 |              6 |                                      0 |                              |                             2 |                         2056 |                         2311 |                         |                         |

## Warnings

- AB: 2 duplicate/missing hour-ending label events in the actual-flow source series; ordered-row alignment is used.
- US: 2 duplicate/missing hour-ending label events in the actual-flow source series; ordered-row alignment is used.

## Claim boundary

Positive published flow is export and negative is import. Source rows are aligned in workbook order because the hour-ending labels contain documented anomalies. TTC is loaded only for the optional hourly-TTC policy; observed-interchange verification has no TTC dependency.
