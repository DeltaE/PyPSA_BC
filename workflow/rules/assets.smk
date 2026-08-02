"""Scenario-independent preparation and validation of generation assets."""


ASSETS = [
    "data/processed_data/wind/existing/bc_ext_wind_assets.csv",
    "data/processed_data/solar/existing/bc_ext_solar_assets.csv",
    "data/processed_data/tpp/existing/bc_ext_tpp_assets.csv",
    "data/processed_data/hydro/existing/hydro_generation.csv",
    "data/processed_data/hydro/existing/hydro_reservoirs.csv",
    "data/processed_data/hydro/existing/hydro_cascade.csv",
    "data/downloaded_data/wind/turbine_dict.json",
]
ASSET_VALIDATION_DIR = f"{OUTPUT_ROOT}/asset_validation"
ASSET_VALIDATION_SUMMARY = f"{ASSET_VALIDATION_DIR}/summary.json"


rule fetch_wind_inventory:
    """Fetch the Canadian wind-turbine inventory used by Rule 2."""
    input:
        config="config/data.yaml",
        script="workflow/scripts/fetch_inputs.py",
    output:
        "data/downloaded_data/wind/canada_turbines.xlsx"
    log:
        "logs/snakemake/02a_fetch_wind_inventory.log"
    run:
        run_logged(
            [PYTHON, "-m", "workflow.scripts.fetch_inputs", "--only", "can_turbines"],
            log[0],
        )


rule asset_preparation:
    """Prepare hydro, wind, solar, and thermal asset inventories."""
    input:
        base=BASE_NETWORK,
        base_validation=BASE_NETWORK_VALIDATION_SUMMARY,
        wind_inventory="data/downloaded_data/wind/canada_turbines.xlsx",
        generators="data/downloaded_data/CODERS/data-pull/supply/generators.csv",
        generic="data/downloaded_data/CODERS/data-pull/supply/generation_generic.csv",
        hydro="data/downloaded_data/CODERS/data-pull/supply/existing_hydro.csv",
        hydro_topology="data/inventory/hydro_topology.csv",
        hydro_reservoirs="data/inventory/hydro_reservoirs.csv",
        generator_wup="data/inventory/custom/hydro_gen_wup_features.csv",
        reservoir_wup="data/inventory/custom/hydro_res_wup_features.csv",
        evidence="studies/chapter2/inputs/cascade_evidence_register.csv",
        # schedules="studies/chapter2/inputs/release_schedules.csv", # to be used in model build
        script="workflow/scripts/prepare_assets.py",
    output:
        ASSETS
    log:
        "logs/snakemake/02_asset_preparation.log"
    run:
        require_validation_pass(input.base_validation, "Base network")
        run_logged([PYTHON, "-m", "workflow.scripts.prepare_assets"], log[0])


rule validate_assets:
    """Retain asset QA evidence and issue a PASS/FAIL result."""
    input:
        assets=ASSETS,
        buses="data/processed_data/network/buses.csv",
        config="config/assets_validation.yaml",
        script="workflow/scripts/validate_assets.py",
    output:
        summary=ASSET_VALIDATION_SUMMARY,
        report=f"{ASSET_VALIDATION_DIR}/validation_report.md",
        schema_issues=f"{ASSET_VALIDATION_DIR}/schema_issues.csv",
        duplicate_identifiers=f"{ASSET_VALIDATION_DIR}/duplicate_identifiers.csv",
        invalid_coordinates=f"{ASSET_VALIDATION_DIR}/invalid_coordinates.csv",
        invalid_capacities=f"{ASSET_VALIDATION_DIR}/invalid_capacities.csv",
        invalid_technologies=f"{ASSET_VALIDATION_DIR}/invalid_technology_labels.csv",
        unmapped_nodes=f"{ASSET_VALIDATION_DIR}/unmapped_network_nodes.csv",
        turbine_coverage=f"{ASSET_VALIDATION_DIR}/wind_turbine_coverage.csv",
        hydro_references=f"{ASSET_VALIDATION_DIR}/hydro_reference_issues.csv",
    log:
        "logs/snakemake/02b_validate_assets.log"
    run:
        run_logged(
            [PYTHON, "-m", "workflow.scripts.validate_assets",
             "--config", input.config,
             "--output-dir", ASSET_VALIDATION_DIR],
            log[0],
        )
