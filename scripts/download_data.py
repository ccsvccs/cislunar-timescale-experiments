#!/usr/bin/env python3
"""Download pinned official DE440 resources; no scientific dependencies needed."""
from pathlib import Path
import argparse
import hashlib
import json
import urllib.request


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def download(name, record, directory, verify_only=False):
    target = directory/name
    if target.is_file():
        if target.stat().st_size != record['bytes'] or digest(target) != record['sha256']:
            raise ValueError(f'Existing file fails pinned size/SHA256: {target}; move it aside before downloading')
        print(f'verified {target}')
        return
    if verify_only:
        raise FileNotFoundError(target)
    part = directory/(name+'.part')
    # 60 seconds is a per-socket timeout, not an unlimited stall; stream in
    # chunks to avoid keeping a large kernel in RAM. Failed partial files are
    # never accepted by the scientific reader.
    request = urllib.request.Request(record['url'], headers={'User-Agent':'tt-tdb-eq21-reproduction/0.1'})
    with urllib.request.urlopen(request, timeout=60) as remote, part.open('wb') as output:
        while block := remote.read(1024*1024):
            output.write(block)
    if part.stat().st_size != record['bytes'] or digest(part) != record['sha256']:
        raise ValueError(f'Download size/SHA256 failed: {part}')
    part.replace(target)
    print(f'downloaded and verified {target}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parents[1]/'data')
    parser.add_argument('--core-only', action='store_true', help='Only de440.bsp and compatible GM PCK')
    parser.add_argument('--verify-only', action='store_true', help='Check files without networking')
    args = parser.parse_args()
    manifest_path = Path(__file__).resolve().parents[1]/'data'/'manifest.json'
    files = json.loads(manifest_path.read_text())['files']
    args.data_dir.mkdir(parents=True, exist_ok=True)
    for name, record in files.items():
        if not args.core_only or record['role'] == 'core':
            download(name, record, args.data_dir, args.verify_only)


if __name__ == '__main__':
    main()
