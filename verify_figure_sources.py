"""Check archived figure inputs, scripts and PDFs without running the analysis."""
from pathlib import Path
import csv
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent


def file_signature(path, mode):
    if mode == 'lf_text':
        content = path.read_bytes().replace(b'\r\n', b'\n')
        return hashlib.sha256(content).hexdigest(), len(content)
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest(), path.stat().st_size


def main():
    manifest = json.loads((ROOT / 'FIGURE_PROVENANCE.json').read_text(encoding='utf-8'))
    records = manifest['raw_inputs'] + manifest['plotting_tables'] + manifest['scripts']
    records += [figure['file'] for figure in manifest['figures']]
    failures = []
    for record in records:
        path = ROOT / record['path']
        if not path.is_file():
            failures.append(f"Missing: {record['path']}")
            continue
        mode = record.get('hash_mode', 'binary')
        checksum, size = file_signature(path, mode)
        expected_size = record.get('canonical_bytes', record['bytes'])
        if size != expected_size or checksum != record['sha256']:
            failures.append(f"Changed: {record['path']}")
        if 'columns' in record:
            with path.open(encoding='utf-8-sig', newline='') as stream:
                columns = next(csv.reader(stream))
            if columns != record['columns']:
                failures.append(f"CSV columns differ: {record['path']}")
    for figure in manifest['figures']:
        for relative in [figure['plotting_script'], *figure['input_tables']]:
            if not (ROOT / relative).is_file():
                failures.append(f"Figure {figure['number']} dependency missing: {relative}")
    if failures:
        print('\n'.join(failures), file=sys.stderr)
        return 1
    print(f"Verified {len(manifest['figures'])} figures, "
          f"{len(manifest['plotting_tables'])} plotting tables, "
          f"{len(manifest['raw_inputs'])} raw inputs and "
          f"{len(manifest['scripts'])} scripts.")
    print('Checks cover archived file integrity and CSV headers; no analysis was run.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
