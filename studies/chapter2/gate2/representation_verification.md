# Chapter 2 hydraulic-representation verification

**Status:** PASS
**Basin-treatment checks passing:** 12/12
**Station-less topology records excluded from transformation:** 1

B and C are generated from the same prepared station-resolved hydro
components. Each row verifies matched storage, annual inflow volume,
turbine water capacity, electrical connection points, explicit external
release routes, removal of obsolete storage-coupling links, and absence
of artificial water self-loops.

| Representation | Cascade | Storage | Inflow | Turbine capacity | Electrical buses | External routes | No self-loops | No placeholders | Status |
|---|---|---|---|---|---|---|---|---|---|
| B_two_store | Bridge | pass | pass | pass | pass | pass | pass | pass | pass |
| B_two_store | Campbell | pass | pass | pass | pass | pass | pass | pass | pass |
| B_two_store | Mica/Columbia | pass | pass | pass | pass | pass | pass | pass | pass |
| B_two_store | Peace | pass | pass | pass | pass | pass | pass | pass | pass |
| B_two_store | Seven Mile | pass | pass | pass | pass | pass | pass | pass | pass |
| B_two_store | Stave | pass | pass | pass | pass | pass | pass | pass | pass |
| C_single_bucket | Bridge | pass | pass | pass | pass | pass | pass | pass | pass |
| C_single_bucket | Campbell | pass | pass | pass | pass | pass | pass | pass | pass |
| C_single_bucket | Mica/Columbia | pass | pass | pass | pass | pass | pass | pass | pass |
| C_single_bucket | Peace | pass | pass | pass | pass | pass | pass | pass | pass |
| C_single_bucket | Seven Mile | pass | pass | pass | pass | pass | pass | pass | pass |
| C_single_bucket | Stave | pass | pass | pass | pass | pass | pass | pass | pass |

These checks establish treatment matching; they do not replace the
unresolved operating-rule evidence gate or full-system validation.
