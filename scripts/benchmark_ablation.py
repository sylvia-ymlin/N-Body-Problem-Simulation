import argparse
import json
import socket

from benchmark_full_suite import ensure_input
from benchmark_suite import run_bench


def run_ablation_suite(binary="./build/nbody_simulate",
                       output_path="data/metrics/uppmax_ablation_results.json"):
    ensure_input("data/inputs/input_50k.gal", 50000, "disk")
    ensure_input("data/inputs/bench_strong_100000.gal", 100000, "disk")

    kmeans_cases = [
        {"label": "morton_no_kmeans", "k": 0},
        {"label": "kmeans_32", "k": 32},
    ]
    schedule_cases = [
        {"label": "static_128", "schedule": "static", "chunk": "128"},
        {"label": "dynamic_128", "schedule": "dynamic", "chunk": "128"},
    ]

    kmeans_results = []
    print("Starting K-Means Ablation...")
    for case in kmeans_cases:
        print(f"  {case['label']}...")
        res = run_bench(
            binary,
            5,
            50000,
            "data/inputs/input_50k.gal",
            1,
            8,
            0.5,
            case["k"],
            extra_env={"NBODY_OMP_SCHEDULE": "dynamic", "NBODY_OMP_CHUNK": "128"},
        )
        if res:
            res.update({"label": case["label"], "k": case["k"], "threads": 8, "N": 50000})
            kmeans_results.append(res)

    schedule_results = []
    print("Starting Schedule Ablation...")
    for case in schedule_cases:
        print(f"  {case['label']}...")
        res = run_bench(
            binary,
            5,
            100000,
            "data/inputs/bench_strong_100000.gal",
            1,
            16,
            0.5,
            0,
            extra_env={
                "NBODY_OMP_SCHEDULE": case["schedule"],
                "NBODY_OMP_CHUNK": case["chunk"],
            },
        )
        if res:
            res.update(
                {
                    "label": case["label"],
                    "threads": 16,
                    "N": 100000,
                    "chunk": int(case["chunk"]),
                }
            )
            schedule_results.append(res)

    payload = {
        "platform": socket.gethostname(),
        "kmeans_ablation": kmeans_results,
        "schedule_ablation": schedule_results,
    }
    with open(output_path, "w") as handle:
        json.dump(payload, handle, indent=4)

    print(f"\nAblation results saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./build/nbody_simulate")
    parser.add_argument("--output", default="data/metrics/uppmax_ablation_results.json")
    args = parser.parse_args()
    run_ablation_suite(binary=args.binary, output_path=args.output)
