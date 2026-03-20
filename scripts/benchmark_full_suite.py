import argparse
import json
import os
import socket
import subprocess

from benchmark_suite import PYTHON_BIN, run_bench, run_comprehensive_suite


def ensure_input(path, n_particles, mode="disk"):
    if os.path.exists(path):
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    subprocess.run(
        [PYTHON_BIN, "scripts/generate_data.py", str(n_particles), path, mode],
        capture_output=True,
        check=True,
    )


def run_diagnostic_suite(binary):
    smoke_file = "data/inputs/smoke_128.gal"
    input_50k = "data/inputs/input_50k.gal"

    ensure_input(smoke_file, 128, "disk")
    ensure_input(input_50k, 50000, "disk")

    cases = [
        {"version": 1, "N": 128, "file": smoke_file, "threads": 1, "theta": 0.0},
        {"version": 2, "N": 128, "file": smoke_file, "threads": 1, "theta": 0.0},
        {"version": 1, "N": 5000, "file": input_50k, "threads": 1, "theta": 0.5},
        {"version": 2, "N": 5000, "file": input_50k, "threads": 1, "theta": 0.5},
        {"version": 1, "N": 10000, "file": input_50k, "threads": 1, "theta": 0.5},
        {"version": 2, "N": 10000, "file": input_50k, "threads": 1, "theta": 0.5},
        {"version": 3, "N": 50000, "file": input_50k, "threads": 1, "theta": 0.5},
        {"version": 4, "N": 50000, "file": input_50k, "threads": 1, "theta": 0.5},
        {"version": 5, "N": 50000, "file": input_50k, "threads": 1, "theta": 0.5},
        {"version": 5, "N": 50000, "file": input_50k, "threads": 16, "theta": 0.5},
    ]

    results = []
    print("Starting Diagnostic Benchmarks...")
    for case in cases:
        print(
            f"  v{case['version']} N={case['N']} threads={case['threads']} theta={case['theta']}..."
        )
        res = run_bench(
            binary,
            case["version"],
            case["N"],
            case["file"],
            1,
            case["threads"],
            case["theta"],
        )
        if res is None:
            continue
        res.update(
            {
                "v": case["version"],
                "n": case["N"],
                "threads": case["threads"],
                "theta": case["theta"],
            }
        )
        results.append(res)
    return results


def run_full_suite(
    binary="./build/nbody_simulate",
    output_path="data/metrics/uppmax_full_results.json",
    threads="1,2,4,8,16",
    n_base=100000,
    theta=0.5,
):
    diagnostic = run_diagnostic_suite(binary)
    scaling_path = output_path + ".scaling.tmp"

    run_comprehensive_suite(
        binary=binary,
        output_path=scaling_path,
        strong_scaling_n=n_base,
        threads_list=[int(part) for part in threads.split(",") if part.strip()],
        nsteps=1,
        theta=theta,
    )

    with open(scaling_path, "r") as handle:
        scaling = json.load(handle)
    os.remove(scaling_path)

    payload = {
        "platform": socket.gethostname(),
        "diagnostic": diagnostic,
        "strong": scaling["strong"],
        "weak": scaling["weak"],
    }

    with open(output_path, "w") as handle:
        json.dump(payload, handle, indent=4)

    print(f"\nFull benchmark results saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./build/nbody_simulate")
    parser.add_argument("--output", default="data/metrics/uppmax_full_results.json")
    parser.add_argument("--threads", default="1,2,4,8,16")
    parser.add_argument("--n-base", type=int, default=100000)
    parser.add_argument("--theta", type=float, default=0.5)
    args = parser.parse_args()

    run_full_suite(
        binary=args.binary,
        output_path=args.output,
        threads=args.threads,
        n_base=args.n_base,
        theta=args.theta,
    )
