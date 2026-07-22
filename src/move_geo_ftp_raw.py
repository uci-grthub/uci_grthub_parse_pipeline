#!/usr/bin/env python3
"""
Relocate the already-uploaded raw FASTQ files on the GEO/NCBI FTP drop box.

The first upload placed the raw FASTQ under FASTQ/<Sublibrary>/<name>. Each GEO
Series needs its raw files flat inside its own experiment folder. NCBI's FTP has
no server-side copy, so this script RENAMES (RNFR/RNTO, no re-transfer) the single
existing copy into the first experiment folder (geo_upload_exp1). The duplicate
copy required by the second Series (geo_upload_exp2) is left to
src/upload_geo_ftp.py, which uploads only the missing files.

After moving, the now-empty FASTQ/<Sublibrary> directories and FASTQ itself are
removed.

Usage:
    python src/move_geo_ftp_raw.py [--dry-run]
"""
import argparse
import ftplib
import logging
import sys
from pathlib import Path

import yaml

PROJECT_DIR = Path(__file__).resolve().parent.parent
FTP_CONFIG_PATH = PROJECT_DIR / "ftp_config.yaml"
CHECKSUMS_FILE = PROJECT_DIR / "data" / "checksums.md5"
DATA_DIR = PROJECT_DIR / "data"
SRC_ROOT = "FASTQ"
DEST_DIR = "geo_upload_exp1"
LOG_FILE = PROJECT_DIR / "logs" / "ftp_move.log"

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


def load_raw_moves():
    """Return (old_rel, new_rel, size) for each raw FASTQ from checksums.md5."""
    moves = []
    with open(CHECKSUMS_FILE) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rel_path = line.split(None, 1)[1]
            local = DATA_DIR / rel_path
            name = Path(rel_path).name
            moves.append((rel_path, f"{DEST_DIR}/{name}", local.stat().st_size))
    return moves


def connect(cfg):
    ftp = ftplib.FTP(cfg["host_address"], timeout=60)
    ftp.login(cfg["username"], cfg["password"])
    ftp.set_pasv(True)
    ftp.cwd("/" + cfg["upload_space"])
    return ftp


def remote_size(ftp, path):
    try:
        ftp.voidcmd("TYPE I")
        return ftp.size(path)
    except ftplib.error_perm:
        return None


def ensure_dest(ftp):
    try:
        ftp.mkd(DEST_DIR)
    except ftplib.error_perm:
        pass  # already exists


def move_one(ftp, old_rel, new_rel, size, dry_run):
    new_size = remote_size(ftp, new_rel)
    if new_size == size:
        log.info("SKIP (already at dest): %s", new_rel)
        return True

    old_size = remote_size(ftp, old_rel)
    if old_size is None:
        log.warning("MISSING source, cannot move: %s", old_rel)
        return False
    if old_size != size:
        log.warning("SIZE MISMATCH source %s: remote=%s local=%s", old_rel, old_size, size)
        return False

    if dry_run:
        log.info("DRY RUN would rename %s -> %s (%d bytes)", old_rel, new_rel, size)
        return True

    ftp.rename(old_rel, new_rel)
    if remote_size(ftp, new_rel) != size:
        log.error("size wrong after rename: %s", new_rel)
        return False
    log.info("MOVED %s -> %s", old_rel, new_rel)
    return True


def remove_empty_src(ftp, dry_run):
    """Remove FASTQ/<Sublibrary> dirs and FASTQ once emptied."""
    try:
        subdirs = ftp.nlst(SRC_ROOT)
    except ftplib.error_perm:
        log.info("no %s tree to clean", SRC_ROOT)
        return
    for entry in subdirs:
        # nlst may return full or relative paths; normalise to path under SRC_ROOT
        sub = entry if entry.startswith(SRC_ROOT) else f"{SRC_ROOT}/{Path(entry).name}"
        if dry_run:
            log.info("DRY RUN would rmd %s", sub)
            continue
        try:
            ftp.rmd(sub)
            log.info("RMD %s", sub)
        except ftplib.error_perm as exc:
            log.warning("could not rmd %s: %s", sub, exc)
    if dry_run:
        log.info("DRY RUN would rmd %s", SRC_ROOT)
        return
    try:
        ftp.rmd(SRC_ROOT)
        log.info("RMD %s", SRC_ROOT)
    except ftplib.error_perm as exc:
        log.warning("could not rmd %s: %s", SRC_ROOT, exc)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="show actions without moving")
    args = parser.parse_args()

    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    cfg = load_ftp_config()
    moves = load_raw_moves()
    log.info("Relocating %d raw files into %s/", len(moves), DEST_DIR)

    ftp = connect(cfg)
    ensure_dest(ftp)
    failures = []
    try:
        for old_rel, new_rel, size in moves:
            if not move_one(ftp, old_rel, new_rel, size, args.dry_run):
                failures.append(old_rel)
        if not failures:
            remove_empty_src(ftp, args.dry_run)
    finally:
        try:
            ftp.quit()
        except Exception:
            pass

    if failures:
        log.error("%d file(s) failed to move:\n%s", len(failures), "\n".join(failures))
        sys.exit(1)
    log.info("All raw files relocated. Run src/upload_geo_ftp.py for the exp2 copy.")


if __name__ == "__main__":
    main()
