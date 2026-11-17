"""Read the Codex canonical JSON; no fallback to parent or DSH mass components."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
_loaded_hashes={}

@lru_cache(maxsize=2)
def _read(name):
    raw=(ROOT/'knowledge'/name).read_bytes()
    current=hashlib.sha256(raw).hexdigest()
    if name in _loaded_hashes and _loaded_hashes[name]!=current:
        raise ValueError('Canonical mass or adapter tables changed in this interpreter; start a fresh process: '+name)
    _loaded_hashes[name]=current
    return json.loads(raw)

def mass_constants():return _read('mass_constants.json')
def components():return _read('components.json')
def formula_text(formula):
    return ''.join(e+(str(n) if n!=1 else '') for e,n in formula.items() if n)
def residue_table():return {n:formula_text(r['residue_formula']) for n,r in components()['residues'].items()}
def delta_table():return {n:(formula_text(r['net_formula']),r['reference_state']) for n,r in components()['covalent_groups'].items()}
def molecule_table():return {n:formula_text(r['formula']) for n,r in components()['neutral_molecules'].items()}
def attachment_table():return {n:set(r['attachment_scope']) for n,r in components()['covalent_groups'].items() if r['scope_status']=='defined_sites'}
def execution_adapters():return _read('execution_adapters.json')

def load_runtime():
    manifest=json.loads((ROOT/'runtime/manifest.json').read_text(encoding='utf-8'))
    for name,expected in manifest['source_sha256'].items():
        actual=hashlib.sha256((ROOT/'knowledge'/name).read_bytes()).hexdigest()
        if actual!=expected:raise ValueError('Stale knowledge runtime; rebuild with tools/build_knowledge.py: '+name)
        if name in _loaded_hashes and _loaded_hashes[name]!=expected:
            raise ValueError('Loaded canonical tables belong to another version; start a fresh process: '+name)
    for name,expected in manifest['compatibility_exports'].items():
        actual=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        if actual!=expected:raise ValueError('Modified derived export; rebuild from canonical knowledge: '+name)
    raw=(ROOT/'runtime/knowledge.json').read_bytes()
    if hashlib.sha256(raw).hexdigest()!=manifest['runtime_sha256']:
        raise ValueError('Runtime checksum mismatch; rebuild tools/build_knowledge.py')
    return json.loads(raw),manifest['logical_sha256']
