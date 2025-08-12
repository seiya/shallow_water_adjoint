import argparse
import os
import resource
import shutil
import subprocess
import time


def run_once(executable: str, threads: int, numa_node: int | None = None) -> tuple[float, float]:
    """Run *executable* with the given number of OpenMP *threads*.

    Optionally bind the process to *numa_node* using ``numactl`` if it is
    available. The function returns a tuple ``(elapsed_seconds, max_rss_mb)``
    where ``max_rss_mb`` is the peak resident set size in MiB.
    """
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["GOMP_CPU_AFFINITY"] = f"0-{threads - 1}"

    cmd = [executable]
    if numa_node is not None and shutil.which("numactl"):
        cpu_range = f"0-{threads - 1}"
        cmd = ["numactl", f"--cpunodebind={numa_node}", f"--physcpubind={cpu_range}", executable]

    start = time.perf_counter()
    subprocess.run(cmd, check=True, env=env, cwd=os.path.dirname(executable))
    elapsed = time.perf_counter() - start

    max_rss_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    max_rss_mb = max_rss_kb / 1024.0
    return elapsed, max_rss_mb


def benchmark(executable: str, thread_list: list[int], numa_node: int | None = None) -> list[tuple[int, float, float]]:
    """Return benchmarking results for *executable*.

    Parameters
    ----------
    executable : str
        Path to the executable to run.
    thread_list : list[int]
        Numbers of threads to benchmark.
    numa_node : int, optional
        If given, bind to the specified NUMA node using ``numactl``.
    """
    results = []
    for t in thread_list:
        elapsed, rss = run_once(executable, t, numa_node)
        results.append((t, elapsed, rss))
    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark an OpenMP executable with optional NUMA binding"
    )
    parser.add_argument("executable", help="Path to executable to benchmark")
    parser.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=[1],
        help="List of OpenMP thread counts to test",
    )
    parser.add_argument(
        "--numa-node",
        type=int,
        default=None,
        help="Bind process to NUMA node using numactl",
    )
    args = parser.parse_args()

    results = benchmark(args.executable, args.threads, args.numa_node)
    print("threads,elapsed_s,max_rss_mb")
    for threads, elapsed, rss in results:
        print(f"{threads},{elapsed:.6f},{rss:.2f}")


if __name__ == "__main__":
    main()
