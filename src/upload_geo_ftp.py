#!/usr/bin/env python3
"""
Upload FASTQ and processed count matrix files to the GEO/NCBI FTP drop box.

Reads connection details from ftp_config.yaml, the raw FASTQ file list from
data/checksums.md5, and the processed files under output/geo_upload_exp1
and output/geo_upload_exp2, then uploads each to the configured upload_space
directory on the FTP server. Skips files whose remote size already matches
the local size, and resumes partial uploads otherwise.

Usage:
    python src/upload_geo_ftp.py [--dry-run] [--retries N]
"""
import argparse
import ftplib
import logging
import sys
import time
from pathlib import Path

import yaml

PROJECT_DIR = Path(__file__).resolve().parent.parent
FTP_CONFIG_PATH = PROJECT_DIR / "ftp_config.yaml"
CHECKSUMS_FILE = PROJECT_DIR / "data" / "checksums.md5"
DATA_DIR = PROJECT_DIR / "data"
PROCESSED_DIRS = ["geo_upload_exp1", "geo_upload_exp2"]
OUTPUT_DIR = PROJECT_DIR / "output"
LOG_FILE = PROJECT_DIR / "logs" / "ftp_upload.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOG_FILE), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)


def load_ftp_config():
    with open(FTP_CONFIG_PATH) as fh:
        cfg = yaml.safe_load(fh)
    for key in ("host_address", "username", "password", "upload_space"):
        if key not in cfg:
            raise ValueError(f"ftp_config.yaml missing required key: {key}")
    return cfg


def load_file_list():
    """Return (local_path, remote_rel_path) pairs for every file to upload.

    Each experiment is a separate GEO Series and must be self-contained, so the
    shared raw FASTQ files are placed (flat, by basename) into every experiment
    folder alongside that experiment's processed files. Remote layout:

        <upload_space>/<experiment>/<basename>

    matching the bare filenames listed in each geo_submission_metadata_*.xlsx.
    """
    raw_paths = []
    with open(CHECKSUMS_FILE) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rel_path = line.split(None, 1)[1]
            raw_paths.append(DATA_DIR / rel_path)

    files = []
    for dirname in PROCESSED_DIRS:
        for raw_local in raw_paths:
            files.append((raw_local, f"{dirname}/{raw_local.name}"))
        for local_path in sorted((OUTPUT_DIR / dirname).iterdir()):
            if not local_path.is_file():
                continue
            if local_path.suffix == ".xlsx" or local_path.name.startswith("sample_list"):
                continue
            files.append((local_path, f"{dirname}/{local_path.name}"))

    return files


def connect(cfg):
    ftp = ftplib.FTP(cfg["host_address"], timeout=60)
    ftp.login(cfg["username"], cfg["password"])
    ftp.set_pasv(True)
    return ftp


def ensure_remote_dir(ftp, remote_dir_parts):
    for part in remote_dir_parts:
        try:
            ftp.cwd(part)
        except ftplib.error_perm:
            ftp.mkd(part)
            ftp.cwd(part)


def remote_size(ftp, filename):
    try:
        ftp.voidcmd("TYPE I")
        return ftp.size(filename)
    except ftplib.error_perm:
        return None


def upload_file(ftp, cfg, local_path, remote_rel, dry_run, retries):
    local_size = local_path.stat().st_size
    parts = Path(remote_rel).parts
    subdirs, filename = parts[:-1], parts[-1]

    for attempt in range(1, retries + 1):
        try:
            ftp.cwd("/" + cfg["upload_space"])
            if subdirs:
                ensure_remote_dir(ftp, subdirs)

            existing = remote_size(ftp, filename)
            if existing == local_size:
                log.info("SKIP (already complete): %s", remote_rel)
                return True

            if dry_run:
                log.info("DRY RUN would upload: %s (%d bytes)", remote_rel, local_size)
                return True

            rest = existing if existing and existing < local_size else None
            with open(local_path, "rb") as fh:
                if rest:
                    fh.seek(rest)
                    log.info("RESUME %s from byte %d", remote_rel, rest)
                else:
                    log.info("UPLOAD %s (%d bytes)", remote_rel, local_size)
                ftp.storbinary(f"STOR {filename}", fh, rest=rest)

            final_size = remote_size(ftp, filename)
            if final_size != local_size:
                raise IOError(
                    f"size mismatch after upload: remote={final_size} local={local_size}"
                )
            log.info("DONE %s", remote_rel)
            return True

        except (*ftplib.all_errors, IOError, OSError) as exc:
            log.warning("attempt %d/%d failed for %s: %s", attempt, retries, remote_rel, exc)
            time.sleep(min(2 ** attempt, 30))
            try:
                ftp.close()
            except Exception:
                pass
            ftp = connect(cfg)

    log.error("FAILED after %d attempts: %s", retries, remote_rel)
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="list actions without uploading")
    parser.add_argument("--retries", type=int, default=5, help="retry attempts per file")
    args = parser.parse_args()

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    cfg = load_ftp_config()
    files = load_file_list()
    log.info("Uploading %d files to %s:%s", len(files), cfg["host_address"], cfg["upload_space"])

    ftp = connect(cfg)
    failures = []
    try:
        for local_path, remote_rel in files:
            ok = upload_file(ftp, cfg, local_path, remote_rel, args.dry_run, args.retries)
            if not ok:
                failures.append(remote_rel)
    finally:
        try:
            ftp.quit()
        except Exception:
            pass

    if failures:
        log.error("%d file(s) failed:\n%s", len(failures), "\n".join(failures))
        sys.exit(1)
    log.info("All files uploaded successfully.")


if __name__ == "__main__":
    main()
