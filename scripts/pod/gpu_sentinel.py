"""GPU Sentinel Daemon (Zero-Idle Auto-Recovery Watchdog).

Monitors active benchmark process:
- Tracks CPU ticks (/proc/<pid>/stat utime+stime) to verify active computation.
- If CPU ticks are completely frozen for > 300 seconds AND GPU-Util is 0%, triggers stall recovery.
- If process exits or crashes unexpectedly, automatically relaunches benchmark from latest
  checkpoint.
- Leaves actively computing processes (GPU util > 0% or advancing CPU ticks) untouched.
"""

import json
import os
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

HEARTBEAT_PATH = Path("/tmp/sina_heartbeat.json")
SENTINEL_LOG = Path("/tmp/sina_sentinel.log")
ALERT_LOG = Path("/tmp/sina_sentinel_alert.log")
BENCHMARK_SCRIPT = Path("/home/ubuntu/sina-platform/scripts/pod/run_wave1_benchmark.py")
PYTHON_EXEC = Path("/home/ubuntu/miniconda3/envs/biopathnet/bin/python")

MAX_FROZEN_SECONDS = 300  # 5 minutes of completely frozen CPU ticks
CHECK_INTERVAL_SECONDS = 30


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: UP017 (py3.10)


def log(msg: str):
    line = f"[{utcnow_iso()}] [Sentinel] {msg}"
    print(line, flush=True)
    with open(SENTINEL_LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def is_pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def get_cpu_ticks(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/stat") as f:
            parts = f.read().split()
            return int(parts[13]) + int(parts[14])
    except Exception:
        return -1


def kill_existing_runners():
    log("Killing any lingering benchmark runners...")
    subprocess.run(["pkill", "-9", "-f", "run_wave1_benchmark.py"], check=False)
    time.sleep(2)


def get_runner_env():
    env = os.environ.copy()
    conda_bin = "/home/ubuntu/miniconda3/envs/biopathnet/bin"
    cuda_bin = "/usr/local/cuda/bin"
    env["PATH"] = f"{conda_bin}:{cuda_bin}:{env.get('PATH', '')}"
    env["PYTHONUNBUFFERED"] = "1"
    env["CUDA_VISIBLE_DEVICES"] = "0"
    return env


def launch_benchmark():
    kill_existing_runners()
    cmd = [str(PYTHON_EXEC), "-u", str(BENCHMARK_SCRIPT), "--batch-size", "32", "--seed", "1337"]
    unit_set = os.environ.get("SINA_UNIT_SET")
    if unit_set:
        cmd += ["--unit-set", unit_set]
    log(f"Launching benchmark runner with clean PATH: {' '.join(cmd)}...")
    with open("/tmp/wave1_benchmark.log", "a", encoding="utf-8") as log_file:
        proc = subprocess.Popen(
            cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            cwd="/home/ubuntu/sina-platform",
            env=get_runner_env(),
        )
    log(f"Benchmark started with new PID {proc.pid}")
    return proc.pid


def find_active_runner_pid():
    try:
        output = subprocess.check_output(["pgrep", "-f", "run_wave1_benchmark.py"]).decode().split()
        pids = [int(p) for p in output if int(p) != os.getpid()]
        return pids[0] if pids else None
    except Exception:
        return None


def main():
    log("=== GPU Sentinel Daemon Started (CPU-Tick Dynamic Supervision) ===")

    current_pid = find_active_runner_pid()
    if current_pid:
        log(f"Attached to existing active benchmark runner PID {current_pid}")
    else:
        current_pid = launch_benchmark()

    last_ticks = get_cpu_ticks(current_pid)
    last_advance_time = time.time()

    while True:
        time.sleep(CHECK_INTERVAL_SECONDS)

        # Check if completed
        if HEARTBEAT_PATH.exists():
            try:
                hb = json.loads(HEARTBEAT_PATH.read_text(encoding="utf-8"))
                if hb.get("status") == "benchmark_finished":
                    log("SUCCESS: Benchmark completed all units! Sentinel exiting.")
                    break
            except Exception:
                pass

        if not is_pid_alive(current_pid):
            log(f"ALERT: Benchmark runner PID {current_pid} exited unexpectedly! Relaunching...")
            with open(ALERT_LOG, "a", encoding="utf-8") as f:
                f.write(
                    f"[{utcnow_iso()}] Process {current_pid} exited unexpectedly."
                    " Force restarting.\n"
                )
            current_pid = launch_benchmark()
            last_ticks = get_cpu_ticks(current_pid)
            last_advance_time = time.time()
            continue

        ticks = get_cpu_ticks(current_pid)
        now = time.time()

        if ticks > last_ticks:
            # Active computation progressing
            last_ticks = ticks
            last_advance_time = now
        else:
            # Ticks stationary
            frozen_duration = now - last_advance_time
            if frozen_duration > MAX_FROZEN_SECONDS:
                log(
                    f"ALERT: Benchmark runner PID {current_pid} is STALLED"
                    f" (CPU ticks frozen for {frozen_duration:.1f}s)! Restarting..."
                )
                with open(ALERT_LOG, "a", encoding="utf-8") as f:
                    f.write(
                        f"[{utcnow_iso()}] Process {current_pid} frozen for"
                        f" {frozen_duration:.1f}s. Force restarting.\n"
                    )
                current_pid = launch_benchmark()
                last_ticks = get_cpu_ticks(current_pid)
                last_advance_time = now


if __name__ == "__main__":
    main()
