"""
Memory monitoring for PL_time benchmark.
Tracks max CPU (RSS) and optionally GPU (PyTorch) memory during execution.
Requires psutil for CPU monitoring.
"""
import os
import threading
import time
from typing import Optional

_max_cpu_gb = [0.0]
_monitor_running = [False]
_psutil_available = [None]


def _check_psutil():
    if _psutil_available[0] is None:
        try:
            import psutil  # noqa: F401
            _psutil_available[0] = True
        except ImportError:
            _psutil_available[0] = False
    return _psutil_available[0]


def _cpu_sample_loop():
    if not _check_psutil():
        return
    import psutil
    p = psutil.Process(os.getpid())
    while _monitor_running[0]:
        try:
            rss = p.memory_info().rss
            _max_cpu_gb[0] = max(_max_cpu_gb[0], rss / (1024**3))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        time.sleep(0.3)


def start_cpu_monitor():
    """Start background thread sampling CPU memory."""
    _max_cpu_gb[0] = 0.0
    _monitor_running[0] = True
    t = threading.Thread(target=_cpu_sample_loop, daemon=True)
    t.start()
    return t


def stop_cpu_monitor(timeout: float = 1.0) -> Optional[float]:
    """Stop monitor and return max CPU memory in GB, or None if psutil unavailable."""
    _monitor_running[0] = False
    time.sleep(0.4)
    if not _check_psutil():
        return None
    return round(_max_cpu_gb[0], 3)


def get_max_gpu_memory_gb() -> Optional[float]:
    """Return max GPU memory allocated (PyTorch) in GB, or None if N/A."""
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.synchronize()
            return round(torch.cuda.max_memory_allocated() / (1024**3), 3)
    except Exception:
        pass
    return None


def reset_gpu_stats():
    """Reset PyTorch peak memory stats (call before GPU work)."""
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
    except Exception:
        pass
