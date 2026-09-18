"""Runnable input -> two stages -> confirmed feedback -> complete export example."""
from pathlib import Path
import sys

import numpy as np

# Allow running this example by absolute path, from any working directory.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lcdl_algorithm import (  # noqa: E402
    Case, SolverOptions, run_two_stage, validate_feedback,
    dispatch_feedback, export_result,
)


def main():
    case = Case.from_json(ROOT / "data/minimal_case.json")
    options = SolverOptions(time_limit=10, max_iterations=5, threads=2)
    result = run_two_stage(case, options)
    print("status:", result["status"])
    print("U:", result["stage1"].get("U"))
    print("L:", result["stage1"].get("L"))

    if result["feasibility_certified"]:
        # Explicit tutorial feedback: first hour 0.44 MW, second hour 0.36 MW.
        # Replace with the user's confirmed K x T full-cycle feedback in an application.
        feedback_p_mw = np.array([[0.44, 0.36]])
        levels = result["stage2"]["levels"]
        check = validate_feedback(case, result["stage1"], levels, feedback_p_mw)
        result["example_feedback_check"] = check
        execution = dispatch_feedback(
            case, result["stage1"], levels, feedback_p_mw, options
        )
        result["example_execution"] = execution
        print("feedback status:", execution["status"])
        if execution["status"] == "optimal":
            dispatch = execution["dispatch"]
            print("charge MW:", dispatch["charge_mw"])
            print("discharge MW:", dispatch["discharge_mw"])
            print("mode z:", dispatch["mode_z"])
            print("energy MWh, including initial:", dispatch["energy_with_initial_mwh"])
            print("cost:", dispatch["cost"])
    else:
        print("Global feasibility not certified; no execution example is accepted.")

    # Also preserves partial results if the global solve has not finished verification.
    output = export_result(
        result, ROOT / "results/readme_example/two_stage_result.json",
        include_matrices=True,
    )
    print("output:", output)


if __name__ == "__main__":
    main()
