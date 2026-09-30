#!/usr/bin/env python3
"""Regenerate the explicit Conda lock and R session record."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-repo",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    args = parser.parse_args()

    benchmark_root = args.benchmark_repo / "benchmarks/brca100_netbid_qc"
    mamba_root = args.benchmark_repo / ".micromamba"
    micromamba = mamba_root / "bin/micromamba"
    environment = mamba_root / "envs/netbid2-r44"
    explicit = subprocess.run(
        [
            str(micromamba), "list", "-r", str(mamba_root),
            "-p", str(environment), "--explicit", "--sha256",
        ],
        text=True,
        capture_output=True,
        check=True,
    ).stdout.splitlines()
    explicit = [line.strip() for line in explicit if line.startswith("http")]
    if not explicit or any(
        re.fullmatch(r"https?://\S+#[0-9a-f]{64}", line) is None
        for line in explicit
    ):
        raise RuntimeError("Micromamba did not return SHA256-pinned package URLs")
    lock_text = "\n".join(
        [
            "# Exact solved Conda packages for linux-64.",
            "# NetBID2 itself is installed and verified by install_netbid2.R.",
            "@EXPLICIT",
            *explicit,
            "",
        ]
    )
    wrapper = benchmark_root / "netbid2-r"
    session_script = benchmark_root / "record_session_info.R"
    session = subprocess.run(
        [str(wrapper), "Rscript", str(session_script)],
        text=True,
        capture_output=True,
        check=True,
    )
    if (
        "NetBID2 version: 2.2.0" not in session.stdout
        or "NetBID2 commit: 5defa454d600b94f5dd6d1f9f4428f99759a6821"
        not in session.stdout
    ):
        raise RuntimeError("Unexpected NetBID2 version/commit in session record")

    lock_path = benchmark_root / "environment-netbid2-r44-linux-64.lock"
    session_path = benchmark_root / "session-info.txt"
    lock_partial = lock_path.with_name(f"{lock_path.name}.{os.getpid()}.partial")
    session_partial = session_path.with_name(
        f"{session_path.name}.{os.getpid()}.partial"
    )
    try:
        lock_partial.write_text(lock_text, encoding="utf-8", newline="\n")
        session_partial.write_text(session.stdout, encoding="utf-8", newline="\n")
        os.replace(lock_partial, lock_path)
        os.replace(session_partial, session_path)
    finally:
        lock_partial.unlink(missing_ok=True)
        session_partial.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
