import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Plot Taylor test convergence")
    parser.add_argument(
        "--test",
        type=int,
        default=1,
        choices=[1, 2],
        help="Taylor test number (1 or 2)",
    )
    args = parser.parse_args()

    build_dir = Path(__file__).resolve().parents[1] / "build"
    data = np.load(build_dir / f"taylor_test{args.test}_results.npz")
    eps = data["eps"]
    fwd = data["fwd"]
    rev = data["rev"]

    plt.loglog(eps, fwd, "o-", label="forward AD")
    plt.loglog(eps, rev, "s-", label="reverse AD")
    plt.xlabel("epsilon")
    plt.ylabel("error")
    plt.legend()
    plt.grid(True, which="both")
    plt.savefig(build_dir / f"taylor_test{args.test}_convergence.png")


if __name__ == "__main__":
    main()
