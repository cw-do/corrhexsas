"""Strict text-file intake. Data stay on the user's computer."""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import numpy as np


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass(frozen=True)
class Data:
    q: object
    intensity: object
    uncertainty: object
    sigma_q: object
    label: str = "data"
    source_hash: str = ""

    def __post_init__(self):
        for name in ("q", "intensity", "uncertainty", "sigma_q"):
            object.__setattr__(self, name, np.asarray(getattr(self, name), float))
        if self.q.ndim != 1 or len(self.q) < 2 or any(getattr(self, k).shape != self.q.shape for k in ("intensity", "uncertainty", "sigma_q")):
            raise ValueError("Data need at least two points with equally sized columns")
        if not np.all(np.isfinite(np.array([self.q, self.intensity, self.uncertainty, self.sigma_q]))):
            raise ValueError("All data columns must be finite")
        if np.any(self.q < 0) or np.any(np.diff(self.q) <= 0):
            raise ValueError("q must be nonnegative and strictly increasing")
        if np.any(self.uncertainty <= 0) or np.any(self.sigma_q < 0):
            raise ValueError("dI must be positive and sigma_q must be nonnegative")


def load_data(path, *, q_unit="1/A", relative_error=None, resolution_relative=0.0,
              resolution_absolute=0.0, label=None):
    """Read q I [dI [sigma_q]] from CSV or whitespace, optional header.

    Two-column files require an explicit relative_error. The configured
    resolution is used only when a fourth column is absent. q_unit applies
    to both the q and sigma_q columns and to resolution_absolute.
    """
    path = Path(path)
    rows, seen_header = [], False
    for number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        words = line.replace(",", " ").split()
        try:
            row = [float(v) for v in words]
        except ValueError as exc:
            if not rows and not seen_header and words[0].lower() in ("q", "q_a-1", "q_a^-1", "q_nm-1"):
                seen_header = True
                continue
            raise ValueError(f"Invalid numeric row in {path.name}, line {number}") from exc
        if len(row) not in (2, 3, 4) or (rows and len(row) != len(rows[0])):
            raise ValueError(f"Expected consistent 2, 3 or 4 columns, line {number}")
        rows.append(row)
    arr = np.asarray(rows, float)
    if arr.ndim != 2 or len(arr) < 2:
        raise ValueError("The file does not contain at least two numeric rows")
    if q_unit not in ("1/A", "1/nm"):
        raise ValueError("q_unit must be 1/A or 1/nm")
    factor = 1.0 if q_unit == "1/A" else 0.1
    q = arr[:, 0]*factor
    if arr.shape[1] >= 3:
        dy = arr[:, 2]
    else:
        if relative_error is None or relative_error <= 0:
            raise ValueError("Two-column input requires an explicit positive relative_error")
        dy = np.abs(arr[:, 1])*relative_error
    if resolution_relative < 0 or resolution_absolute < 0:
        raise ValueError("Resolution widths must be nonnegative")
    dq = arr[:, 3]*factor if arr.shape[1] == 4 else np.hypot(resolution_relative*q, resolution_absolute*factor)
    return Data(q, arr[:, 1], dy, dq, label or path.stem, sha256(path))
