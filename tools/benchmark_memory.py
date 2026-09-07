#!/usr/bin/env python3
"""Measure peak resident memory for an SRFM command and its child processes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import time

import psutil


def measure_command(command, interval=0.05):
    """Execute a command while sampling aggregate resident memory.

    Args:
        command (list[str]): Command and arguments to execute.
        interval (float): Sampling interval in seconds.

    Returns:
        dict: Exit status, runtime, peak RSS, and periodic RSS samples.
    """
    started = time.perf_counter()
    process = subprocess.Popen(command)
    monitored = psutil.Process(process.pid)
    peak_rss = 0
    peak_root_rss = 0
    peak_child_rss = 0
    samples = []
    while process.poll() is None:
        resident = 0
        root_resident = 0
        child_resident = 0
        try:
            try:
                root_resident = monitored.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            for child in monitored.children(recursive=True):
                try:
                    child_resident += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        except psutil.NoSuchProcess:
            pass
        resident = root_resident + child_resident
        elapsed = time.perf_counter() - started
        peak_rss = max(peak_rss, resident)
        peak_root_rss = max(peak_root_rss, root_resident)
        peak_child_rss = max(peak_child_rss, child_resident)
        samples.append(
            {
                "seconds": elapsed,
                "root_rss_bytes": root_resident,
                "child_rss_bytes": child_resident,
                "aggregate_rss_bytes": resident,
            }
        )
        time.sleep(interval)
    return {
        "command": command,
        "returncode": process.returncode,
        "elapsed_seconds": time.perf_counter() - started,
        "peak_rss_bytes": peak_rss,
        "peak_root_rss_bytes": peak_root_rss,
        "peak_child_rss_bytes": peak_child_rss,
        "samples": samples,
    }


def main():
    """Parse command-line arguments and write a JSON memory report.

    Returns:
        None: The process exits with the measured command's status.
    """
    parser = argparse.ArgumentParser(
        description="Sample aggregate RSS for an SRFM command and its subprocesses."
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval", type=float, default=0.05)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    arguments = parser.parse_args()
    command = arguments.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command must follow --")
    report = measure_command(command, arguments.interval)
    arguments.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    raise SystemExit(report["returncode"])


if __name__ == "__main__":
    main()
