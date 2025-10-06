import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional


def _find_results_path(test_number: int, user_dir: Optional[Path]) -> Path:
    filename = f"taylor_test{test_number}_results.npz"

    if user_dir is not None:
        path = (user_dir / filename).resolve()
        if not path.exists():
            raise FileNotFoundError(
                f"Results file '{filename}' not found in provided directory '{user_dir}'."
            )
        return path

    repo_root = Path(__file__).resolve().parents[1]
    candidate_dirs = [
        Path.cwd(),
        repo_root / "nompi" / "build",
        repo_root / "mpi" / "build",
    ]

    for directory in candidate_dirs:
        path = (directory / filename).resolve()
        if path.exists():
            return path

    searched = ", ".join(str(directory) for directory in candidate_dirs)
    raise FileNotFoundError(
        f"Results file '{filename}' not found. Searched directories: {searched}. "
        "Pass --results-dir to specify a custom location."
    )


def main():
    parser = argparse.ArgumentParser(description="Plot Taylor test convergence")
    parser.add_argument(
        "--test",
        type=int,
        default=1,
        choices=[1, 2, 5],
        help="Taylor test number (1, 2, or 5)",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        help="Directory containing taylor test results (defaults to common build directories)",
    )
    args = parser.parse_args()

    results_path = _find_results_path(args.test, args.results_dir)
    data = np.load(results_path)
    eps = data["eps"]
    fwd = data["fwd"]
    rev = data["rev"]

    plt.loglog(eps, fwd, "o-", label="forward AD")
    plt.loglog(eps, rev, "s-", label="reverse AD")
    plt.xlabel("epsilon")
    plt.ylabel("error")
    plt.legend()
    plt.grid(True, which="both")
    output_path = results_path.with_name(f"taylor_test{args.test}_convergence.png")
    plt.savefig(output_path)


if __name__ == "__main__":
    main()
