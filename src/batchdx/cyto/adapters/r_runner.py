"""Running R from Python without embedding an R interpreter in this process.

Every R call is a subprocess running a script from ``workflow/scripts``.  Inputs and
outputs cross the boundary as JSON plus the typed binary matrix format in
:mod:`batchdx.cyto.io.binary`, so an R-side failure is a non-zero exit status with a
captured log rather than a corrupted in-process state.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from batchdx.cyto.schemas.contracts import SoftwareVersion


class RNotAvailable(RuntimeError):
    """Raised when no usable R interpreter is on PATH."""


class RScriptFailed(RuntimeError):
    """Raised when an R script exits non-zero.  Carries the captured log."""

    def __init__(self, script: str, returncode: int, log: str):
        super().__init__(f"{script} exited with status {returncode}")
        self.script = script
        self.returncode = returncode
        self.log = log


@dataclass(frozen=True)
class RResult:
    """Outcome of one R subprocess call."""

    returncode: int
    stdout: str
    stderr: str
    payload: dict

    @property
    def log(self) -> str:
        return self.stdout + ("\n" + self.stderr if self.stderr else "")


def r_available(rscript: str = "Rscript") -> bool:
    return shutil.which(rscript) is not None


def r_version(rscript: str = "Rscript") -> SoftwareVersion:
    """Resolve the R version actually on PATH, or report it as unavailable."""

    if not r_available(rscript):
        return SoftwareVersion(
            name="R", available=False, source="PATH",
            detail=f"{rscript!r} was not found on PATH",
        )
    proc = subprocess.run(
        [rscript, "-e", "cat(paste(R.version$major, R.version$minor, sep='.'))"],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        return SoftwareVersion(name="R", available=False, source="PATH",
                               detail=proc.stderr.strip()[:500])
    return SoftwareVersion(name="R", version=proc.stdout.strip(), available=True, source="PATH")


def r_package_version(package: str, rscript: str = "Rscript") -> SoftwareVersion:
    """Report an R package's installed version without assuming it is present."""

    if not r_available(rscript):
        return SoftwareVersion(name=package, available=False, source="R library",
                               detail="no R interpreter on PATH")
    code = (
        f"v <- tryCatch(as.character(utils::packageVersion({package!r})), "
        f"error=function(e) NA_character_); cat(ifelse(is.na(v), '', v))"
    )
    proc = subprocess.run([rscript, "-e", code], capture_output=True, text=True, check=False)
    version = proc.stdout.strip()
    if proc.returncode != 0 or not version:
        return SoftwareVersion(
            name=package, available=False, source="R library",
            detail=(proc.stderr.strip() or f"package {package!r} is not installed")[:500],
        )
    return SoftwareVersion(name=package, version=version, available=True, source="R library")


def script_path(name: str) -> Path:
    """Locate a packaged R script, whether running from a checkout or an installed wheel."""

    from batchdx.cyto.paths import workflow_root

    candidate = workflow_root() / "scripts" / name
    if not candidate.exists():
        raise FileNotFoundError(f"packaged R script not found: {candidate}")
    return candidate


def run_r_script(
    name: str,
    params: dict,
    work_dir: str | Path,
    *,
    rscript: str = "Rscript",
    log_path: str | Path | None = None,
    timeout: float | None = None,
) -> RResult:
    """Run ``workflow/scripts/<name>`` with ``params`` serialised to JSON.

    The script writes its structured result to ``<work_dir>/result.json``; anything it
    prints is captured into the log.  A non-zero exit raises :class:`RScriptFailed`,
    never a silently unchanged matrix.
    """

    if not r_available(rscript):
        raise RNotAvailable(f"{rscript!r} is not on PATH")

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    params_file = work / "params.json"
    result_file = work / "result.json"
    if result_file.exists():
        result_file.unlink()
    params_file.write_text(json.dumps(params, indent=2, default=str) + "\n", encoding="utf-8")

    script = script_path(name)
    proc = subprocess.run(
        [rscript, "--vanilla", str(script), str(params_file), str(result_file)],
        capture_output=True, text=True, check=False, timeout=timeout,
    )
    log = proc.stdout + ("\n" + proc.stderr if proc.stderr else "")
    if log_path is not None:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        Path(log_path).write_text(log, encoding="utf-8")

    if proc.returncode != 0:
        raise RScriptFailed(name, proc.returncode, log)

    payload: dict = {}
    if result_file.exists():
        payload = json.loads(result_file.read_text(encoding="utf-8"))
    return RResult(returncode=proc.returncode, stdout=proc.stdout,
                   stderr=proc.stderr, payload=payload)
