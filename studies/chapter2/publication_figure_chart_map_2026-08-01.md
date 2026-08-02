# Chapter 2 publication figure chart map

This source note defines the analytical contract for current publication figures.
It is not manuscript prose.

| Figure | Analytical question | Evidence and grain | Chart contract | Encoding and scale | Output and QA surface | Publication role |
|---|---|---|---|---|---|---|
| Named-reservoir turbine cost and local marginal price | Does economic ordering explain why named turbines are bypassed? | Six named reservoir links; calendar-2021 effective VOM, 8,760-hour mean and maximum local marginal price | Comparison; horizontal paired dot/lollipop | Zero-based CAD/MWh axis; dark filled circle for mean price, blue open circle for maximum, orange square for VOM; shape distinguishes series without color | `results/workflow/sensitivities/hydro_dispatch_economics/figure_hydro_dispatch_economics.png` and `.svg`; inspected at exported resolution | Mechanism diagnostic; not observed-price validation |
| Named-station hydro dispatch sensitivity | How does imposing station-energy bands redistribute hydro generation and water routing? | Seven named stations; two calendar-2021 model solutions plus one fiscal-2021 comparator | Two-panel comparison: paired dot plot for GWh and grouped horizontal bars for turbine capture | Both axes start at zero; reference uses dark circle/bar, sensitivity blue square/bar, fiscal comparator orange open circle; markers and fill states support grayscale reading | `results/workflow/sensitivities/hydro_dispatch_comparison/figure_hydro_dispatch_sensitivity.png` and `.svg`; inspected at exported resolution | Primary hydro-policy sensitivity figure |
| Annual generation mix calibration | Does the reference reproduce declared provincial generation categories? | Annual calendar-2021 model and Statistics Canada Table 25-10-0015-01 | Category comparison | TWh units and explicit evidence-role notes | `results/workflow/calibration/annual_reference_2021/figure_generation_mix.png` and `.svg` | Calibration fidelity; not held-out validation |
| Monthly hydro calibration | Does modeled hydro reproduce annual magnitude and monthly shape? | Twelve calendar months | Highlighted two-series line | Common TWh/month axis; 12 temporal points | `results/workflow/calibration/annual_reference_2021/figure_monthly_hydro.png` and `.svg` | Calibration fidelity |
| Load-input fidelity | Does regional disaggregation preserve the source provincial chronology? | 8,760 hourly observations plus duration/summary comparisons | Trend and benchmark comparison | GW units; source-versus-model distinction; source chronology is an input | `results/workflow/calibration/annual_reference_2021/figure_load_input_fidelity.png` and `.svg` | Input-fidelity check only |
| Annual dispatch | Is the solved reference internally balanced without emergency backstop? | 8,760 hourly dispatch records | Stacked composition over time plus balance diagnostics | GW units; positive supply components; backstop explicit | `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_diagnostics/figure_annual_dispatch.png` and `.svg` | Computational integrity |
| Reservoir storage | Do stores remain within bounds and satisfy cyclic closure? | Hourly storage trajectories by reservoir | Small-multiple/selected-series trend | Fraction or volume units must stay explicit; bounds visible | `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_diagnostics/figure_reservoir_storage.png` and `.svg` | Computational integrity and mechanism context |
| Corridor utilization | Which modeled transport corridors carry stress after circulation removal? | Hourly lexicographic minimum-transfer reconstruction | Ranked comparison/duration diagnostic | Per-unit utilization from zero; proxy-capacity limitation adjacent to figure | `results/workflow/annual/A_lower_zero_observed_transport_full_horizon_diagnostics/figure_corridor_utilization.png` and `.svg` | Diagnostic only; no validated congestion claim |

## Current visual QA result

The two new hydro figures pass the static visual review: titles and units are
visible, axes are honest, labels and legends do not clip, the palette is restrained,
and marker shape or fill supplements color. Their captions state the period and
evidence boundary. The existing calibration and integrity figures remain eligible
for the supplement, but each must be inspected in the final manuscript or HTML
container before submission.

## Missing primary figures

The publication package still lacks empirical A/B/C representation effects because
the readiness gate remains on hold. The following figures must not be generated from
unapproved scenario runs:

1. paired A/B/C effect plot across matched scenario families;
2. stress-event chronology linking water state, releases, congestion, and shortfall;
3. routing-delay sensitivity;
4. release-rule and terminal-storage sensitivity.

The input network and cascade schematics already exist under `vis/input_visuals/`,
but their final manuscript selection requires a separate source, resolution, and
caption audit.
