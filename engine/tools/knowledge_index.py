"""Read-only reference queries, with explicit stale-index detection."""
import json
from pathlib import Path
import sqlite3
import hashlib
from contextlib import contextmanager
from knowledge_source import ROOT, load_runtime
@contextmanager
def open_index():
    _,version=load_runtime()
    p=ROOT/'index/knowledge.sqlite3'
    manifest=json.loads((ROOT/'index/manifest.json').read_text(encoding='utf-8'))
    if hashlib.sha256(p.read_bytes()).hexdigest()!=manifest['database_sha256']:
        raise ValueError('Index checksum mismatch; rebuild tools/build_knowledge_index.py')
    con=sqlite3.connect(p.as_uri()+'?mode=ro',uri=True)
    source=con.execute("SELECT value FROM metadata WHERE key='source_logical_sha256'").fetchone()
    if not source or source[0]!=version:
        con.close();raise ValueError('Stale index; run tools/build_knowledge_index.py')
    try:yield con
    finally:con.close()
def lookup(rule_id):
    with open_index() as con:
        row=con.execute('SELECT record_json FROM rules WHERE rule_id=?',(rule_id,)).fetchone()
        return json.loads(row[0]) if row else None
def mass_leads(value,tolerance):
    if tolerance<0:raise ValueError('Mass tolerance must be nonnegative')
    with open_index() as con:
        return [json.loads(r[0]) for r in con.execute('SELECT record_json FROM rules WHERE delta_mono BETWEEN ? AND ? ORDER BY rule_id',(value-tolerance,value+tolerance))]
def templates(kind=None,value=None):
    with open_index() as con:
        if kind is None:return [json.loads(r[0]) for r in con.execute('SELECT recipe_json FROM branches ORDER BY branch_id')]
        return [json.loads(r[0]) for r in con.execute('SELECT b.recipe_json FROM branches b JOIN memberships m USING(branch_id) WHERE m.kind=? AND m.value=? ORDER BY b.branch_id',(kind,value))]
def find_records(text):
    with open_index() as con:
        return [json.loads(r[0]) for r in con.execute('SELECT record_json FROM rules WHERE record_json LIKE ? ORDER BY rule_id',('%'+text+'%',))]
