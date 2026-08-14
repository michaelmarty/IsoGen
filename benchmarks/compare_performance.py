"""Compare MSVC and Intel wheel benchmark results."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msvc", type=Path, required=True)
    parser.add_argument("--intel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--maximum-slowdown", type=float, default=1.10)
    arguments = parser.parse_args()

    msvc = json.loads(arguments.msvc.read_text(encoding="utf-8"))
    intel = json.loads(arguments.intel.read_text(encoding="utf-8"))
    if msvc["benchmarks"].keys() != intel["benchmarks"].keys():
        raise RuntimeError("The benchmark files contain different test cases")

    comparisons = {}
    ratios = []
    lines = [
        "## Windows compiler performance",
        "",
        "| Workload | MSVC (us) | Intel (us) | Intel / MSVC |",
        "| --- | ---: | ---: | ---: |",
    ]
    for name in msvc["benchmarks"]:
        msvc_case = msvc["benchmarks"][name]
        intel_case = intel["benchmarks"][name]
        np.testing.assert_allclose(
            intel_case["result"],
            msvc_case["result"],
            rtol=2e-5,
            atol=2e-6,
            err_msg=f"Numerical results differ for {name}",
        )
        msvc_seconds = msvc_case["median_seconds"]
        intel_seconds = intel_case["median_seconds"]
        ratio = intel_seconds / msvc_seconds
        ratios.append(ratio)
        comparisons[name] = {
            "msvc_seconds": msvc_seconds,
            "intel_seconds": intel_seconds,
            "intel_to_msvc_ratio": ratio,
        }
        lines.append(
            f"| {name} | {msvc_seconds * 1e6:.2f} | "
            f"{intel_seconds * 1e6:.2f} | {ratio:.3f} |"
        )

    geometric_mean = math.exp(sum(math.log(ratio) for ratio in ratios) / len(ratios))
    passed = geometric_mean <= arguments.maximum_slowdown
    lines.extend(
        [
            "",
            f"Geometric-mean Intel/MSVC ratio: **{geometric_mean:.3f}**.",
            f"Allowed maximum: **{arguments.maximum_slowdown:.3f}**.",
        ]
    )
    report = {
        "geometric_mean_intel_to_msvc_ratio": geometric_mean,
        "maximum_slowdown": arguments.maximum_slowdown,
        "passed": passed,
        "comparisons": comparisons,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    summary = "\n".join(lines) + "\n"
    print(summary)
    if arguments.summary:
        with arguments.summary.open("a", encoding="utf-8") as stream:
            stream.write(summary)

    if not passed:
        print("Intel wheel exceeded the allowed aggregate slowdown")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
