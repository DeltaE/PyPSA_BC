"""Scenario-independent preparation and validation of the electrical base network."""


BASE_NETWORK = [
    "data/processed_data/network/buses.csv",
    "data/processed_data/network/lines.csv",
    "data/processed_data/network/line_types.csv",
    "data/processed_data/network/transformers.csv",
    "data/processed_data/network/transformer_types.csv",
]
BASE_NETWORK_VALIDATION_DIR = f"{OUTPUT_ROOT}/base_network_validation"
BASE_NETWORK_VALIDATION_SUMMARY = f"{BASE_NETWORK_VALIDATION_DIR}/summary.json"


rule base_network_preparation:
    """Create detailed substation-level buses and transmission tables."""
    input:
        data_config="config/data.yaml",
        assumptions="config/base_network.yaml",
        lines="data/downloaded_data/CODERS/data-pull/network/transmission_lines.csv",
        substations="data/downloaded_data/CODERS/data-pull/network/substations.csv",
        generators="data/downloaded_data/CODERS/data-pull/supply/generators.csv",
        line_table="data/inventory/custom/electric_power_generation_table_13_3a.xlsx",
        script="workflow/scripts/prepare_base_network.py",
    output:
        BASE_NETWORK
    log:
        "logs/snakemake/01_base_network_preparation.log"
    run:
        run_logged([PYTHON, "-m", "workflow.scripts.prepare_base_network"], log[0])


rule validate_base_network:
    """Preserve topology and parameter evidence and issue a PASS/FAIL result."""
    input:
        network=BASE_NETWORK,
        config="config/base_network_validation.yaml",
        script="workflow/scripts/validate_base_network.py",
    output:
        summary=f"{BASE_NETWORK_VALIDATION_DIR}/summary.json",
        report=f"{BASE_NETWORK_VALIDATION_DIR}/validation_report.md",
        invalid_buses=f"{BASE_NETWORK_VALIDATION_DIR}/invalid_buses.csv",
        missing_endpoints=f"{BASE_NETWORK_VALIDATION_DIR}/missing_line_endpoints.csv",
        self_loops=f"{BASE_NETWORK_VALIDATION_DIR}/self_loops.csv",
        components=f"{BASE_NETWORK_VALIDATION_DIR}/disconnected_components.csv",
        parameter_outliers=f"{BASE_NETWORK_VALIDATION_DIR}/parameter_outliers.csv",
        duplicate_identifiers=f"{BASE_NETWORK_VALIDATION_DIR}/duplicate_identifiers.csv",
        parallel_lines=f"{BASE_NETWORK_VALIDATION_DIR}/parallel_lines.csv",
    log:
        "logs/snakemake/01b_validate_base_network.log"
    run:
        run_logged(
            [PYTHON, "-m", "workflow.scripts.validate_base_network",
             "--network-dir", "data/processed_data/network",
             "--config", input.config,
             "--output-dir", BASE_NETWORK_VALIDATION_DIR],
            log[0],
        )
