"""Shared run-directory, provenance, download, and audit utilities."""

from __future__ import annotations

import contextlib
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_SUBDIRS = ("raw", "predictions", "derived", "figures", "logs", "audit")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def directory_size(path: Path) -> int:
    return sum(
        item.stat().st_size
        for item in path.rglob("*")
        if item.is_file() and not item.name.startswith("._")
    )


def peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    # macOS reports bytes; Linux and most BSD environments report KiB.
    return value if sys.platform == "darwin" else value * 1024


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), *args], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def checkout_revision(path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def discover_alphagenome_source() -> Path | None:
    """Find the SDK checkout used by the development workspace, if present."""
    for parent in (REPO_ROOT, *REPO_ROOT.parents):
        candidate = parent / "alphagenome" / "src"
        if candidate.is_dir():
            return candidate
    return None


def make_alphagenome_importable() -> Path | None:
    try:
        __import__("alphagenome")
        return None
    except ModuleNotFoundError:
        source = discover_alphagenome_source()
        if source is not None and str(source) not in sys.path:
            sys.path.insert(0, str(source))
        return source


def load_env_file(path: Path) -> list[str]:
    """Load a small KEY=VALUE file without printing or returning secret values."""
    loaded: list[str] = []
    if not path.exists():
        return loaded
    for line_number, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[7:].strip()
        if "=" not in line:
            raise ValueError(f"Malformed environment line {line_number} in {path}")
        name, value = line.split("=", 1)
        name = name.strip()
        value = value.strip()
        if value[:1] == value[-1:] and value[:1] in {"'", '"'}:
            value = value[1:-1]
        if not name.replace("_", "").isalnum() or name[:1].isdigit():
            raise ValueError(f"Invalid environment name on line {line_number} in {path}")
        if name not in os.environ:
            os.environ[name] = value
            loaded.append(name)
    return loaded


def load_api_key_file(path: Path | None) -> bool:
    """Load a raw key or one-line ALPHAGENOME_API_KEY assignment without logging it."""
    if path is None:
        return False
    value = path.read_text().strip()
    if not value or "\n" in value or "\r" in value:
        raise ValueError(f"API key file must contain exactly one non-empty line: {path}")
    if "=" in value:
        name, value = value.split("=", 1)
        if name.strip() != "ALPHAGENOME_API_KEY":
            raise ValueError(
                f"Key-file assignment must use ALPHAGENOME_API_KEY: {path}"
            )
        value = value.strip()
        if value[:1] == value[-1:] and value[:1] in {"'", '"'}:
            value = value[1:-1]
        if not value:
            raise ValueError(f"ALPHAGENOME_API_KEY is empty in {path}")
    if "ALPHAGENOME_API_KEY" not in os.environ:
        os.environ["ALPHAGENOME_API_KEY"] = value
        return True
    return False


def api_key() -> str:
    value = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not credential_present():
        raise RuntimeError(
            "Missing or placeholder ALPHAGENOME_API_KEY. Copy .env.example to .env and add the key, "
            "or export it in the shell."
        )
    return value


def credential_present() -> bool:
    value = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    lowered = value.lower()
    return bool(value) and "replace_with" not in lowered and "your_authorized_key" not in lowered


def initialize_run(run_dir: Path, *, resume: bool) -> None:
    run_dir = run_dir.resolve()
    if run_dir.exists() and any(run_dir.iterdir()) and not resume:
        raise FileExistsError(
            f"Run directory is not empty: {run_dir}. Choose a new directory or pass --resume."
        )
    run_dir.mkdir(parents=True, exist_ok=True)
    for name in RUN_SUBDIRS:
        (run_dir / name).mkdir(exist_ok=True)


def archive_previous_audit(run_dir: Path) -> None:
    """Retain the previous attempt before a resumed command updates audit files."""
    audit_dir = run_dir / "audit"
    current = audit_dir / "run.json"
    if not current.exists():
        return
    try:
        previous = json.loads(current.read_text())
        stamp = str(previous.get("started_utc", utc_now())).replace(":", "").replace("+", "_")
    except (OSError, json.JSONDecodeError):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    destination = audit_dir / "attempts" / stamp
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(current, destination / "run.json")
    report = audit_dir / "REPRODUCIBILITY_REPORT.md"
    if report.exists():
        shutil.copy2(report, destination / report.name)
    comparison = audit_dir / "comparison.json"
    if comparison.exists():
        shutil.copy2(comparison, destination / comparison.name)
        comparison.unlink()


def download(url: str, destination: Path, *, resume: bool = True) -> dict[str, object]:
    """Download atomically and return non-secret provenance."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and resume:
        return {
            "url": url,
            "path": str(destination),
            "bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
            "reused": True,
        }
    partial = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "SORT1-reanalysis/1.0"})
    last_error: Exception | None = None
    attempts = 3
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            partial.replace(destination)
            break
        except Exception as exc:
            last_error = exc
            if attempt < attempts:
                time.sleep(2 * attempt)
    else:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Download failed after {attempts} attempts: {url}: {last_error}") from last_error
    return {
        "url": url,
        "path": str(destination),
        "bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
        "reused": False,
        "attempts": attempt,
    }


@dataclass
class Audit:
    run_dir: Path
    panels: list[str]
    started_utc: str = field(default_factory=utc_now)
    steps: list[dict[str, object]] = field(default_factory=list)
    downloads: list[dict[str, object]] = field(default_factory=list)
    api_calls: dict[str, int] = field(default_factory=dict)
    api_requests: dict[str, int] = field(default_factory=dict)
    status: str = "running"
    error: str | None = None

    @property
    def json_path(self) -> Path:
        return self.run_dir / "audit" / "run.json"

    def base_record(self) -> dict[str, object]:
        sdk_source = discover_alphagenome_source()
        return {
            "schema_version": 1,
            "status": self.status,
            "error": self.error,
            "started_utc": self.started_utc,
            "updated_utc": utc_now(),
            "panels": self.panels,
            "run_directory": str(self.run_dir.resolve()),
            "repository": {
                "path": str(REPO_ROOT),
                "commit": git_value("rev-parse", "HEAD"),
                "branch": git_value("branch", "--show-current"),
                "dirty": bool(git_value("status", "--porcelain")),
            },
            "environment": {
                "python": sys.version.split()[0],
                "platform": platform.platform(),
                "packages": {
                    name: package_version(name)
                    for name in (
                        "alphagenome",
                        "alphagenome-research",
                        "numpy",
                        "pandas",
                        "pysam",
                        "pyliftover",
                        "matplotlib",
                    )
                },
                "api_key_present": credential_present(),
                "alphagenome_source_checkout": str(sdk_source) if sdk_source else None,
                "alphagenome_source_commit": checkout_revision(
                    sdk_source.parent if sdk_source else None
                ),
            },
            "command": sys.argv,
            "reproduction_code": reproduction_code_manifest(),
            "steps": self.steps,
            "downloads": self.downloads,
            "api_calls": self.api_calls,
            "api_requests": self.api_requests,
            "disk_bytes": directory_size(self.run_dir),
            "peak_rss_bytes": peak_rss_bytes(),
        }

    def save(self) -> None:
        self.json_path.parent.mkdir(parents=True, exist_ok=True)
        self.json_path.write_text(json.dumps(self.base_record(), indent=2) + "\n")

    @contextlib.contextmanager
    def step(self, name: str) -> Iterator[dict[str, object]]:
        record: dict[str, object] = {"name": name, "started_utc": utc_now(), "status": "running"}
        self.steps.append(record)
        self.save()
        start = time.perf_counter()
        try:
            yield record
        except BaseException as exc:
            record.update(
                status="failed",
                elapsed_seconds=round(time.perf_counter() - start, 3),
                error=f"{type(exc).__name__}: {exc}",
            )
            self.status = "failed"
            self.error = str(record["error"])
            self.save()
            raise
        else:
            record.update(status="complete", elapsed_seconds=round(time.perf_counter() - start, 3))
            self.save()

    def add_api_calls(self, panel: str, count: int) -> None:
        """Record scored prediction units (variants or complete sequences)."""
        self.api_calls[panel] = self.api_calls.get(panel, 0) + int(count)
        self.save()

    def add_api_requests(self, panel: str, count: int = 1) -> None:
        """Record actual client method invocations separately from prediction units."""
        self.api_requests[panel] = self.api_requests.get(panel, 0) + int(count)
        self.save()

    def finish(self) -> None:
        self.status = "complete"
        self.error = None
        self.save()


def file_manifest(root: Path) -> list[dict[str, object]]:
    rows = []
    for path in sorted(root.rglob("*")):
        if (
            path.is_file()
            and not path.name.startswith("._")
            and path.name not in {"run.json", "comparison.json", "REPRODUCIBILITY_REPORT.md"}
        ):
            rows.append(
                {
                    "path": str(path.relative_to(root)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    return rows


def reproduction_code_manifest() -> list[dict[str, object]]:
    paths = [REPO_ROOT / "reproduce.py", *sorted((REPO_ROOT / "reproduction").glob("*.py"))]
    return [
        {
            "path": str(path.relative_to(REPO_ROOT)),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in paths
        if path.is_file()
    ]
