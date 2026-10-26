"""Versioned built-in XGBoost models and installation identities."""
from pathlib import Path
import hashlib
import json
import shutil
import tarfile
import tempfile
import urllib.request


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def builtin_models():
    root = Path(__file__).resolve().parents[1] / 'resources/xgboost'
    return {record['name']: record for path in sorted(root.glob('*.json'))
        for record in [json.loads(path.read_text())]}


def model_identity(model_f, name):
    model_f = Path(model_f)
    checksum = sha256(model_f)
    manifest_f = model_f.with_name(model_f.name.removesuffix('_xgboost_model.json') + '_manifest.json')
    model_id = name
    if manifest_f.exists():
        manifest = json.loads(manifest_f.read_text())
        expected = manifest['files'][model_f.name]
        if checksum != expected:
            raise ValueError(f'XGBoost model checksum mismatch: {model_f}')
        model_id = manifest['model_id']
    return {'model_id': model_id, 'model_sha256': checksum}


def install_builtin(name, dest_dir):
    manifest = builtin_models()[name]
    dest_dir = Path(dest_dir)
    dest_dir.parent.mkdir(parents=True, exist_ok=True)
    # Validate the complete download before replacing any installed files.
    with tempfile.TemporaryDirectory(dir=dest_dir.parent) as directory:
        stage = Path(directory)
        archive = stage / 'download.tar.gz'
        urllib.request.urlretrieve(manifest['url'], archive)
        if sha256(archive) != manifest['archive_sha256']:
            raise ValueError(f'Archive checksum mismatch for {name}')
        found = set()
        with tarfile.open(archive, 'r:gz') as tar:
            for member in tar.getmembers():
                path = Path(member.name)
                if member.isdir() and path.parts == (manifest['model_id'],):
                    continue
                if (not member.isfile() or len(path.parts) != 2 or
                        path.parts[0] != manifest['model_id'] or
                        path.name not in manifest['files'] or path.name in found):
                    raise ValueError(f'Unexpected model archive member: {member.name}')
                with tar.extractfile(member) as source, (stage / path.name).open('wb') as target:
                    shutil.copyfileobj(source, target)
                if sha256(stage / path.name) != manifest['files'][path.name]:
                    raise ValueError(f'Model member checksum mismatch: {path.name}')
                found.add(path.name)
        if found != set(manifest['files']):
            raise ValueError(f'Incomplete model archive for {name}')
        dest_dir.mkdir(exist_ok=True)
        for filename in found:
            shutil.copy2(stage / filename, dest_dir / filename)
        (dest_dir / f'{name}_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


def resolve_xgboost_entry(entry, scdata_dir):
    import polars as pl
    registry = pl.read_csv(Path(scdata_dir) / 'xgboost/xgboost_models.csv')
    rows = registry.filter(pl.col('model') == entry['model'])
    if len(rows) != 1:
        raise ValueError(f"Expected one registered XGBoost model named {entry['model']}; found {len(rows)}")
    row = rows.row(0, named=True)
    for key in ['model_f', 'cls_f', 'genes_f']:
        entry[key] = row[key]
        if not Path(entry[key]).is_file():
            raise FileNotFoundError(f"XGBoost file missing: {entry[key]}; rerun scprocess setup")
    entry['label_map_f'] = row.get('label_map_f') or ''
    entry.update(model_identity(entry['model_f'], entry['model']))
    return entry
