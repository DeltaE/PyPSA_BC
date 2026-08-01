"""Shared execution and validation helpers for the PyPSA-BC workflow."""

import json
import os
import subprocess


def run_logged(arguments, logfile):
    """Run a command in the repository, mirror output, and retain a complete log."""
    log_path = Path(str(logfile))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    path_parts = [str(REPO / "src"), str(REPO)]
    if env.get("PYTHONPATH"):
        path_parts.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(path_parts)
    completed = subprocess.run(
        [str(item) for item in arguments],
        cwd=REPO,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(completed.stdout or "", encoding="utf-8")
    if completed.stdout:
        print(completed.stdout, end="")
    if completed.returncode:
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}; see {log_path}"
        )


def run_configured_scenario(
    case_file,
    network_file,
    report_file,
    gate1_file,
    gate2_file,
    representation_gate_file,
    release_uncertainty_gate_file,
    logfile,
):
    """Apply the solve authorization and invoke the scenario adapter."""
    if not config["execution"]["allow_solve"]:
        raise RuntimeError(
            "Production solve is disabled in config/workflow.yaml. Pass Gate 1, review "
            "the case file, then set execution.allow_solve: true."
        )
    for gate_path in (
        gate1_file,
        gate2_file,
        representation_gate_file,
        release_uncertainty_gate_file,
    ):
        gate = json.loads(Path(str(gate_path)).read_text(encoding="utf-8"))
        if gate.get("status") != "PASS":
            raise RuntimeError(
                f"Scientific gate is not PASS: {gate_path} "
                f"(status={gate.get('status', 'missing')})."
            )
    arguments = [
        PYTHON, "-m", "workflow.scripts.run_scenario",
        "--scenario", case_file,
        "--network", network_file,
        "--report", report_file,
        "--year", str(config["execution"]["year"]),
        "--capacity-choice", config["execution"]["capacity_choice"],
    ]
    if config["execution"]["include_vre_investments"]:
        arguments.append("--include-vre-investments")
    run_logged(arguments, logfile)


def run_base_model_solve(prepared_network, solved_network, summary, gates, logfile):
    """Enforce solve authorization and all scientific gates before optimization."""
    if not config["execution"]["allow_solve"]:
        raise RuntimeError(
            "Base-model solve is disabled. Inspect the prepared network and validation "
            "outputs, then set execution.allow_solve: true in config/workflow.yaml."
        )
    for gate_path in gates:
        gate = json.loads(Path(str(gate_path)).read_text(encoding="utf-8"))
        if gate.get("status") != "PASS":
            raise RuntimeError(
                f"Base-model solve blocked by {gate_path}: "
                f"status={gate.get('status', 'missing')}."
            )
    run_logged(
        [PYTHON, "-m", "workflow.scripts.solve_model",
         "--network", prepared_network, "--output", solved_network,
         "--summary", summary],
        logfile,
    )


def run_base_model_preparation(network, report, representation, water_policy,
                               release_multiplier, logfile):
    """Invoke the build-only adapter with the reviewed base-model settings."""
    arguments = [
        PYTHON, "-m", "workflow.scripts.prepare_model_file",
        "--network", network,
        "--report", report,
        "--year", str(config["execution"]["year"]),
        "--capacity-choice", config["execution"]["capacity_choice"],
        "--reservoir-representation", representation,
        "--water-policy", water_policy,
        "--release-multiplier", str(release_multiplier),
    ]
    if config["execution"]["include_vre_investments"]:
        arguments.append("--include-vre-investments")
    run_logged(arguments, logfile)


def require_validation_pass(summary_file, label):
    """Stop downstream work while retaining the validator's evidence artifacts."""
    summary = json.loads(Path(str(summary_file)).read_text(encoding="utf-8"))
    if summary.get("status") != "PASS":
        blockers = ", ".join(summary.get("blocking_checks", [])) or "unspecified"
        raise RuntimeError(
            f"{label} validation is not PASS; blocking checks: {blockers}. "
            f"Inspect {summary_file} and its companion CSV files."
        )
