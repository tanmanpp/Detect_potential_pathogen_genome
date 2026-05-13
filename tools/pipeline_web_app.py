#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Local web API for Detect_potential_pathogen_pipeline.py.

The browser submits a FASTQ file and pipeline options to this server. The server
then runs the existing command-line pipeline either natively or through WSL.
"""

from __future__ import annotations

import json
import mimetypes
import os
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
import uuid
from dataclasses import dataclass, field
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path, PureWindowsPath
from urllib.parse import quote, unquote, urlparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = PROJECT_ROOT / "tools"
WEB_DIR = PROJECT_ROOT / "web"
RUNS_DIR = PROJECT_ROOT / "web_runs"
JOBS_DIR = RUNS_DIR / "jobs"
RESULTS_DIR = RUNS_DIR / "results"
DEFAULT_RCF_ENV = os.environ.get("PATHOGEN_RCF_ENV", "DPPP-rcf")

PIPELINE_SCRIPT = TOOLS_DIR / "Detect_potential_pathogen_pipeline.py"
SUMMARY_SCRIPT = TOOLS_DIR / "generate_pipeline_summary_html.py"
KRAKEN_HTML_SCRIPT = TOOLS_DIR / "kraken2_html_excel_report.py"
MAPPING_SCRIPT = TOOLS_DIR / "mapping_multi_genomes.py"
BLAST_SCRIPT = TOOLS_DIR / "run_blastn_with_virulent_database.sh"

DEFAULT_RUNTIME = os.environ.get("PATHOGEN_WEB_RUNTIME", "auto").lower()
DEFAULT_WSL_DISTRO = os.environ.get("PATHOGEN_WEB_WSL_DISTRO", "")


@dataclass
class Job:
    id: str
    sample_id: str
    created_at: float
    job_dir: Path
    outdir: Path
    command: list[str]
    runtime: str
    process: subprocess.Popen | None = None
    status: str = "queued"
    returncode: int | None = None
    error: str | None = None
    lock: threading.Lock = field(default_factory=threading.Lock)

    @property
    def log_path(self) -> Path:
        return self.job_dir / "web_runner.log"

    @property
    def sample_dir(self) -> Path:
        return self.outdir / self.sample_id

    @property
    def final_dir(self) -> Path:
        return self.sample_dir / "final_results"

    @property
    def pipeline_log_path(self) -> Path:
        return self.sample_dir / "logs" / f"{self.sample_id}.pipeline.log"


JOBS: dict[str, Job] = {}
JOBS_LOCK = threading.Lock()


def sanitize_id(value: str) -> str:
    value = value.strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("._-") or "sample"


def choose_runtime(requested: str) -> str:
    requested = (requested or DEFAULT_RUNTIME).lower()
    if requested == "auto":
        return "wsl" if os.name == "nt" and shutil.which("wsl") else "native"
    if requested in {"wsl", "native"}:
        return requested
    raise ValueError("runtime must be auto, wsl, or native")


def windows_path_to_wsl(path: Path | str) -> str:
    raw = str(path)
    if raw.startswith("/"):
        return raw
    win = PureWindowsPath(raw)
    drive = win.drive.rstrip(":").lower()
    if not drive:
        return raw.replace("\\", "/")
    rest = "/".join(win.parts[1:])
    return f"/mnt/{drive}/{rest}"


def wsl_path_to_windows(path: Path | str) -> str:
    raw = str(path).strip()
    match = re.match(r"^/mnt/([A-Za-z])/(.*)$", raw)
    if not match:
        return raw
    drive = match.group(1).upper()
    rest = match.group(2).replace("/", "\\")
    return f"{drive}:\\{rest}" if rest else f"{drive}:\\"


def looks_like_windows_path(value: str) -> bool:
    return bool(re.match(r"^[A-Za-z]:[\\/]", value.strip()))


def clean_user_path(value: str) -> str:
    value = value.strip()
    while len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        value = value[1:-1].strip()
    return value


def runtime_path(value: Path | str, runtime: str) -> str:
    text = str(value)
    if runtime == "wsl":
        if looks_like_windows_path(text) or (os.name == "nt" and not text.startswith("/")):
            return windows_path_to_wsl(text)
    return text


def user_path(value: str, runtime: str) -> str:
    value = clean_user_path(value)
    if not value:
        return value
    if looks_like_windows_path(value):
        if runtime == "wsl" or os.name != "nt":
            return windows_path_to_wsl(value)
        return value
    if value.startswith("/mnt/") and runtime == "native" and os.name == "nt":
        return wsl_path_to_windows(value)
    return value


def resolve_output_dir(value: str, runtime: str) -> Path:
    value = clean_user_path(value)
    if not value:
        return RESULTS_DIR.resolve()

    # On WSL/Linux hosts, accept Windows-style absolute paths from the UI.
    if os.name != "nt" and looks_like_windows_path(value):
        return Path(windows_path_to_wsl(value)).expanduser().resolve()

    if runtime == "wsl":
        # If the UI sends a Windows path while the pipeline runs in WSL,
        # keep it as a local Windows path on Windows hosts, but convert it to
        # a WSL mount path on POSIX hosts.
        if looks_like_windows_path(value):
            if os.name == "nt":
                return Path(value).expanduser().resolve()
            return Path(windows_path_to_wsl(value)).expanduser().resolve()

        # If the server itself runs on Windows and the UI sends a WSL mount
        # path, convert it back so the local server can create directories.
        if os.name == "nt" and value.startswith("/mnt/"):
            return Path(wsl_path_to_windows(value)).expanduser().resolve()

    return Path(value).expanduser().resolve()


def parse_bool(value: str | None) -> bool:
    return str(value or "").lower() in {"1", "true", "yes", "on"}


def get_field(fields: dict[str, str], name: str, default: str = "") -> str:
    value = fields.get(name, default)
    return value.strip() if isinstance(value, str) else default


def build_pipeline_args(fields: dict[str, str], reads_path: Path, outdir: Path, runtime: str) -> list[str]:
    sample_id = sanitize_id(get_field(fields, "sample_id"))
    human_index = clean_user_path(get_field(fields, "human_index"))
    if not human_index:
        raise ValueError("human_index is required")

    threads = get_field(fields, "threads", "8")
    min_len = get_field(fields, "min_len", "200")
    kraken_sc = get_field(fields, "kraken_sc", "0.2")

    args = [
        runtime_path(PIPELINE_SCRIPT, runtime),
        "-i", runtime_path(reads_path, runtime),
        "-sid", sample_id,
        "-o", runtime_path(outdir, runtime),
        "-t", threads,
        "--min_len", min_len,
        "--human_index", user_path(human_index, runtime),
        "--kraken_db", user_path(get_field(fields, "kraken_db", "/home/yilun/YiLun/kraken_database/standard"), runtime),
        "--rcf_env", get_field(fields, "rcf_env", DEFAULT_RCF_ENV),
        "--rcf_taxdb", user_path(get_field(fields, "rcf_taxdb", "/home/yilun/YiLun/kraken_database/standard"), runtime),
        "--kraken_sc", kraken_sc,
        "--kraken_html_script", runtime_path(KRAKEN_HTML_SCRIPT, runtime),
        "--mapping_script", runtime_path(MAPPING_SCRIPT, runtime),
        "--vg_bacteria", runtime_path(PROJECT_ROOT / "Bacteria", runtime),
        "--vg_virus", runtime_path(PROJECT_ROOT / "Virus", runtime),
        "--vg_parasite", runtime_path(PROJECT_ROOT / "Parasite", runtime),
        "--blast_script", runtime_path(BLAST_SCRIPT, runtime),
        "--vfdb", runtime_path(PROJECT_ROOT / "VFDB" / "VFDB_setA_nt_rename.fas", runtime),
    ]

    for flag in [
        "skip_len_filter",
        "skip_rm_human",
        "skip_flye",
        "skip_kraken",
        "skip_mapping",
        "skip_blast",
    ]:
        if parse_bool(fields.get(flag)):
            args.append(f"--{flag}")

    return args


def build_command(pipeline_args: list[str], runtime: str) -> list[str]:
    if runtime == "wsl":
        shell_cmd = shlex.join(["python3", *pipeline_args])
        cmd = ["wsl"]
        if DEFAULT_WSL_DISTRO:
            cmd.extend(["-d", DEFAULT_WSL_DISTRO])
        cmd.extend(["bash", "-lc", shell_cmd])
        return cmd
    return [sys.executable, *pipeline_args]


def write_job_state(job: Job) -> None:
    payload = {
        "id": job.id,
        "sample_id": job.sample_id,
        "created_at": job.created_at,
        "status": job.status,
        "returncode": job.returncode,
        "runtime": job.runtime,
        "outdir": str(job.outdir),
        "command": job.command,
        "error": job.error,
    }
    job.job_dir.mkdir(parents=True, exist_ok=True)
    (job.job_dir / "job.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def start_job(job: Job) -> None:
    def runner() -> None:
        with job.lock:
            job.status = "running"
            write_job_state(job)

        with job.log_path.open("a", encoding="utf-8", errors="replace") as log:
            log.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting job {job.id}\n")
            log.write(f"Runtime: {job.runtime}\n")
            log.write("Command:\n")
            log.write(shlex.join(job.command) + "\n\n")
            log.flush()
            try:
                job.process = subprocess.Popen(
                    job.command,
                    cwd=str(PROJECT_ROOT),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                rc = job.process.wait()
                with job.lock:
                    job.returncode = rc
                    job.status = "finished" if rc == 0 else "failed"
            except Exception as exc:
                with job.lock:
                    job.error = str(exc)
                    job.status = "failed"
                log.write(f"\n[ERROR] {exc}\n")

        write_job_state(job)

    threading.Thread(target=runner, daemon=True).start()


def parse_multipart(body: bytes, content_type: str) -> tuple[dict[str, str], dict[str, tuple[str, bytes]]]:
    match = re.search(r"boundary=(?P<boundary>[^;]+)", content_type)
    if not match:
        raise ValueError("Missing multipart boundary")
    boundary = match.group("boundary").strip('"').encode()
    fields: dict[str, str] = {}
    files: dict[str, tuple[str, bytes]] = {}

    for part in body.split(b"--" + boundary):
        part = part.strip()
        if not part or part == b"--":
            continue
        if part.endswith(b"--"):
            part = part[:-2].strip()
        header_blob, _, data = part.partition(b"\r\n\r\n")
        if not header_blob:
            continue
        data = data.rstrip(b"\r\n")
        headers = header_blob.decode("utf-8", errors="replace").split("\r\n")
        disposition = next((h for h in headers if h.lower().startswith("content-disposition:")), "")
        name_match = re.search(r'name="([^"]+)"', disposition)
        if not name_match:
            continue
        name = name_match.group(1)
        filename_match = re.search(r'filename="([^"]*)"', disposition)
        if filename_match:
            filename = Path(filename_match.group(1)).name or "reads.fastq"
            files[name] = (filename, data)
        else:
            fields[name] = data.decode("utf-8", errors="replace")

    return fields, files


def tail_file(path: Path, max_bytes: int = 120_000) -> str:
    if not path.exists():
        return ""
    size = path.stat().st_size
    with path.open("rb") as fh:
        if size > max_bytes:
            fh.seek(size - max_bytes)
        data = fh.read()
    return data.decode("utf-8", errors="replace")


def list_result_files(job: Job) -> list[dict[str, str | int]]:
    if not job.final_dir.exists():
        return []
    out = []
    for fp in sorted(job.final_dir.iterdir()):
        if fp.is_file():
            out.append({
                "name": fp.name,
                "size": fp.stat().st_size,
                "url": f"/api/jobs/{quote(job.id)}/files/{quote(fp.name)}",
            })
    return out


def job_payload(job: Job) -> dict:
    if job.process is not None and job.status == "running":
        rc = job.process.poll()
        if rc is not None:
            job.returncode = rc
            job.status = "finished" if rc == 0 else "failed"
            write_job_state(job)

    summary_files = [f for f in list_result_files(job) if str(f["name"]).endswith(".summary_report.html")]
    return {
        "id": job.id,
        "sample_id": job.sample_id,
        "status": job.status,
        "returncode": job.returncode,
        "runtime": job.runtime,
        "created_at": job.created_at,
        "outdir": str(job.outdir),
        "sample_dir": str(job.sample_dir),
        "final_dir": str(job.final_dir),
        "error": job.error,
        "logs": {
            "web_runner": tail_file(job.log_path),
            "pipeline": tail_file(job.pipeline_log_path),
        },
        "files": list_result_files(job),
        "summary_url": summary_files[0]["url"] if summary_files else None,
    }


class Handler(BaseHTTPRequestHandler):
    server_version = "PathogenPipelineWeb/0.1"

    def log_message(self, fmt: str, *args) -> None:
        sys.stderr.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), fmt % args))

    def send_json(self, payload: dict, status: int = 200) -> None:
        data = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def send_file(self, path: Path, download_name: str | None = None) -> None:
        if not path.exists() or not path.is_file():
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        ctype = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        data = path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        if download_name:
            self.send_header("Content-Disposition", f'inline; filename="{download_name}"')
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        path = parsed.path

        if path in {"/", "/index.html"}:
            self.send_file(WEB_DIR / "index.html")
            return

        if path == "/api/health":
            self.send_json({
                "ok": True,
                "project_root": str(PROJECT_ROOT),
                "runtime_default": DEFAULT_RUNTIME,
                "wsl_available": bool(shutil.which("wsl")),
            })
            return

        if path == "/api/jobs":
            with JOBS_LOCK:
                jobs = [job_payload(job) for job in JOBS.values()]
            jobs.sort(key=lambda x: x["created_at"], reverse=True)
            self.send_json({"jobs": jobs})
            return

        match = re.fullmatch(r"/api/jobs/([^/]+)", path)
        if match:
            job = JOBS.get(unquote(match.group(1)))
            if not job:
                self.send_error(HTTPStatus.NOT_FOUND)
                return
            self.send_json(job_payload(job))
            return

        match = re.fullmatch(r"/api/jobs/([^/]+)/files/([^/]+)", path)
        if match:
            job = JOBS.get(unquote(match.group(1)))
            if not job:
                self.send_error(HTTPStatus.NOT_FOUND)
                return
            filename = Path(unquote(match.group(2))).name
            self.send_file(job.final_dir / filename, filename)
            return

        self.send_error(HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)

        if parsed.path == "/api/jobs":
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if length <= 0:
                    raise ValueError("Empty request")
                fields, files = parse_multipart(self.rfile.read(length), self.headers.get("Content-Type", ""))
                if "reads" not in files:
                    raise ValueError("Please choose a FASTQ/FASTA reads file")

                runtime = choose_runtime(get_field(fields, "runtime", DEFAULT_RUNTIME))
                sample_id = sanitize_id(get_field(fields, "sample_id"))
                if not sample_id:
                    raise ValueError("sample_id is required")

                job_id = f"{time.strftime('%Y%m%d-%H%M%S')}-{sample_id}-{uuid.uuid4().hex[:8]}"
                job_dir = JOBS_DIR / job_id
                upload_dir = job_dir / "uploads"
                upload_dir.mkdir(parents=True, exist_ok=True)
                filename, data = files["reads"]
                reads_path = upload_dir / Path(filename).name
                reads_path.write_bytes(data)

                outdir_text = get_field(fields, "outdir")
                outdir = resolve_output_dir(outdir_text, runtime)
                outdir.mkdir(parents=True, exist_ok=True)

                pipeline_args = build_pipeline_args(fields, reads_path, outdir, runtime)
                command = build_command(pipeline_args, runtime)

                job = Job(
                    id=job_id,
                    sample_id=sample_id,
                    created_at=time.time(),
                    job_dir=job_dir,
                    outdir=outdir,
                    command=command,
                    runtime=runtime,
                )
                with JOBS_LOCK:
                    JOBS[job.id] = job
                write_job_state(job)
                start_job(job)
                self.send_json(job_payload(job), status=201)
            except Exception as exc:
                self.send_json({"error": str(exc)}, status=400)
            return

        match = re.fullmatch(r"/api/jobs/([^/]+)/cancel", parsed.path)
        if match:
            job = JOBS.get(unquote(match.group(1)))
            if not job:
                self.send_error(HTTPStatus.NOT_FOUND)
                return
            if job.process and job.status == "running":
                job.process.terminate()
                job.status = "cancelled"
                write_job_state(job)
            self.send_json(job_payload(job))
            return

        self.send_error(HTTPStatus.NOT_FOUND)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Run the local web UI for the pathogen detection pipeline.")
    parser.add_argument("--host", default="127.0.0.1", help="Host/IP to bind")
    parser.add_argument("--port", type=int, default=8080, help="Port to bind")
    args = parser.parse_args()

    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    url = f"http://{args.host}:{args.port}"
    print(f"[OK] Pathogen pipeline web UI is running: {url}")
    print("[INFO] Keep this terminal open while jobs are running.")

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[INFO] Server stopped.")


if __name__ == "__main__":
    main()
