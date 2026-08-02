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
BASE_NETWORK_CORRECTION_AUDIT = "data/processed_data/network/correction_audit.csv"


rule fetch_base_network_validation:
    """Fetch and fingerprint the public map used by registered manual corrections."""
    input:
        config="config/data.yaml",
        script="workflow/scripts/fetch_base_network_validation.py",
    output:
        map="data/validation/base_network/bc_hydro_transmission_system_2025.pdf",
        metadata="data/validation/base_network/bc_hydro_transmission_system_2025.metadata.json",
    log:
        "logs/snakemake/01a_fetch_base_network_validation.log"
    run:
        run_logged(
            [PYTHON, "-m", "workflow.scripts.fetch_base_network_validation"],
            log[0],
        )


rule base_network_preparation:
    """Create detailed substation-level buses and transmission tables."""
    input:
        data_config="config/data.yaml",
        assumptions="config/base_network.yaml",
        lines="data/downloaded_data/CODERS/data-pull/network/transmission_lines.csv",
        substations="data/downloaded_data/CODERS/data-pull/network/substations.csv",
        generators="data/downloaded_data/CODERS/data-pull/supply/generators.csv",
        line_table="data/inventory/custom/electric_power_generation_table_13_3a.xlsx",
        node_corrections="data/validation/base_network/node_corrections.csv",
        line_corrections="data/validation/base_network/line_corrections.csv",
        map="data/validation/base_network/bc_hydro_transmission_system_2025.pdf",
        map_metadata="data/validation/base_network/bc_hydro_transmission_system_2025.metadata.json",
        script="workflow/scripts/prepare_base_network.py",
    output:
        BASE_NETWORK + [BASE_NETWORK_CORRECTION_AUDIT]
    log:
        "logs/snakemake/01_base_network_preparation.log"
    run:
        run_logged([PYTHON, "-m", "workflow.scripts.prepare_base_network"], log[0])


rule validate_base_network:
    """Preserve topology and parameter evidence and issue a PASS/FAIL result."""
    input:
        network=BASE_NETWORK,
        correction_audit=BASE_NETWORK_CORRECTION_AUDIT,
        config="config/base_network_validation.yaml",
        component_policy="data/validation/base_network/component_policy.yaml",
        node_corrections="data/validation/base_network/node_corrections.csv",
        line_corrections="data/validation/base_network/line_corrections.csv",
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
