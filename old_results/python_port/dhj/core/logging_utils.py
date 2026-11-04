"""Tee stdout/stderr into a run.log file, line-buffered, so the log on disk
always reflects what has been printed so far (a killed job keeps its log)."""
import datetime as _dt
import os
import sys


class _Tee:
    def __init__(self, stream, fh):
        self._stream = stream
        self._fh = fh

    def write(self, s):
        self._stream.write(s)
        self._fh.write(s)
        if "\n" in s:
            self._fh.flush()

    def flush(self):
        self._stream.flush()
        self._fh.flush()

    def isatty(self):
        return False

    def fileno(self):
        return self._stream.fileno()

    @property
    def encoding(self):
        return getattr(self._stream, "encoding", "utf-8")


class RunLogger:
    """Context manager that mirrors everything printed to <out_dir>/run.log.

    The file is opened in append mode so a resumed run continues the same log;
    a banner with the timestamp and the command line separates the sessions.
    """

    def __init__(self, out_dir: str, filename: str = "run.log"):
        os.makedirs(out_dir, exist_ok=True)
        self.path = os.path.join(out_dir, filename)
        self._fh = None
        self._saved = None

    def __enter__(self):
        self._fh = open(self.path, "a", encoding="utf-8", buffering=1)
        self._saved = (sys.stdout, sys.stderr)
        sys.stdout = _Tee(self._saved[0], self._fh)
        sys.stderr = _Tee(self._saved[1], self._fh)
        stamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("#" * 70)
        print(f"# run started {stamp}")
        print(f"# command: {' '.join(sys.argv)}")
        print(f"# cwd: {os.getcwd()}")
        print("#" * 70)
        return self

    def __exit__(self, exc_type, exc, tb):
        stamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if exc_type is not None:
            print(f"# run FAILED {stamp}: {exc_type.__name__}: {exc}")
        else:
            print(f"# run finished {stamp}")
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout, sys.stderr = self._saved
        self._fh.close()
        return False
