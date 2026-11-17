"""Codex-native deterministic SQLite index. Dynamic edits are indexed as templates, not fake deltas."""
import hashlib
import json
from pathlib import Path
import re
import sqlite3
from knowledge_source import ROOT, load_runtime

def encoded(x):return json.dumps(x,sort_keys=True,ensure_ascii=False,separators=(',',':'))
def digest(x):return hashlib.sha256(encoded(x).encode()).hexdigest()
def build():
    bundle,source_digest=load_runtime();constants=bundle['mass_constants']
    def mono(f):return sum(constants['monoisotopic_atoms'][e]*n for e,n in f.items())
    def average(f):return sum(constants['average_atoms'][e]*n for e,n in f.items())
    rule_rows=[];branch_rows=[];memberships=[];component_rows=[];appendix_rows=[]
    for r in bundle['rules']:
        legacy=r['legacy_record'] or {};f=legacy.get('delta_formula')
        static=mono(f) if isinstance(f,dict) else 0.0 if legacy.get('mass_status')=='isomeric' else None
        nominal=re.match(r'^\s*([+−–-]?\d+)\s*amu',str(legacy.get('nominal_mass_shift_book','')))
        n=int(nominal[1].replace('−','-').replace('–','-')) if nominal else None
        d=r['disposition'] or {}
        rule_rows.append((r['rule_id'],r['origin'],d.get('status','legacy_reference_lead'),static,
                          encoded(f) if f is not None else None,n,encoded(r),average(f) if isinstance(f,dict) else 0.0 if legacy.get('mass_status')=='isomeric' else None))
        for b in r['branches']:
            recipe=b['recipe'];transform=recipe['transform'];selector=recipe.get('site_selector',recipe.get('bond_selector',{}))
            # Only explicitly defined fixed deltas are mass-indexed. Other branches instantiate on state.
            delta=transform.get('delta_formula')
            value=mono(delta) if isinstance(delta,dict) else None
            branch_rows.append((b['branch_id'],r['rule_id'],transform['kind'],value,
                                encoded(delta) if delta is not None else None,encoded(recipe),average(delta) if isinstance(delta,dict) else None))
            values={'residue':set(selector.get('residues',[]))|set(selector.get('requires_residues',[]))|set(selector.get('motif',[])),
                'condition':set(recipe.get('required_conditions',[]))|set(recipe.get('required_any_conditions',[]))|set(recipe.get('supporting_any_conditions',[])),
                'group':{selector[k] for k in ('group','required_group') if isinstance(selector.get(k),str)},
                'bond':set(recipe.get('bond_selector',{}).get('kinds',[])),
                'source_figure':set(str(x) for x in recipe['source_refs'].get('figures',[]))}
            for category,items in values.items():
                for item in sorted(items):memberships.append((category,item,b['branch_id']))
    comp=bundle['components']
    for name,r in comp['residues'].items():
        for kind,key in [('residue','residue_formula'),('free_amino_acid','free_amino_acid_formula')]:
            f=r[key];component_rows.append((kind+':'+name,kind,name,mono(f),encoded(f),encoded(r),average(f)))
    for cat,key in [('covalent_groups','net_formula'),('neutral_molecules','formula')]:
        for name,r in comp[cat].items():
            f=r[key];component_rows.append((cat+':'+name,cat,name,mono(f),encoded(f),encoded(r),average(f)))
    for row in bundle['appendix_dispositions']:
        appendix_rows.append((row['rule_id'],row['status'],row['book_raw_value'],row['book_mass_basis'],
                              mono(row['delta_formula']),row['intrinsic_charge'],encoded(row),average(row['delta_formula'])))
    tables={'rules':sorted(rule_rows),'branches':sorted(branch_rows),'memberships':sorted(set(memberships)),
            'components':sorted(component_rows),'appendix':sorted(appendix_rows)}
    logical=digest(tables)
    dest=ROOT/'index/knowledge.sqlite3';dest.parent.mkdir(exist_ok=True)
    temp=dest.with_suffix('.building.sqlite3')
    assert ROOT in temp.resolve().parents
    if temp.exists():temp.unlink()
    con=sqlite3.connect(temp);con.execute('PRAGMA foreign_keys=ON')
    con.executescript('''
    CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
    CREATE TABLE rules(rule_id TEXT PRIMARY KEY,origin TEXT,status TEXT,delta_mono REAL,delta_formula TEXT,book_nominal INTEGER,record_json TEXT NOT NULL,delta_average REAL);
    CREATE TABLE branches(branch_id TEXT PRIMARY KEY,rule_id TEXT NOT NULL REFERENCES rules(rule_id),transform_kind TEXT,delta_mono REAL,delta_formula TEXT,recipe_json TEXT NOT NULL,delta_average REAL);
    CREATE TABLE memberships(kind TEXT,value TEXT,branch_id TEXT REFERENCES branches(branch_id),PRIMARY KEY(kind,value,branch_id));
    CREATE TABLE components(component_id TEXT PRIMARY KEY,kind TEXT,name TEXT,mono REAL,formula TEXT,record_json TEXT,average_mass REAL);
    CREATE TABLE appendix(row_id TEXT PRIMARY KEY,status TEXT,book_value INTEGER,book_basis TEXT,element_delta_mono REAL,intrinsic_charge INTEGER,record_json TEXT,element_delta_average REAL);
    CREATE INDEX rules_mass ON rules(delta_mono); CREATE INDEX rules_nominal ON rules(book_nominal);
    CREATE INDEX branches_mass ON branches(delta_mono); CREATE INDEX branches_transform ON branches(transform_kind);
    CREATE INDEX component_mass ON components(mono); CREATE INDEX appendix_book ON appendix(book_value);
    CREATE INDEX rules_average ON rules(delta_average); CREATE INDEX branches_average ON branches(delta_average);
    ''')
    con.executemany('INSERT INTO metadata VALUES (?,?)',sorted({'schema':'codex-index-1','source_logical_sha256':source_digest,'index_logical_sha256':logical,
                                                              'compiler_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}.items()))
    for name,rows in tables.items():
        if rows:con.executemany(f"INSERT INTO {name} VALUES ({','.join('?' for _ in rows[0])})",rows)
    assert not con.execute('PRAGMA foreign_key_check').fetchall()
    checks=0
    for table,rows,col_idx,id_idx,mass_column in [('rules',tables['rules'],3,0,'delta_mono'),('branches',tables['branches'],3,0,'delta_mono'),
                                               ('rules',tables['rules'],7,0,'delta_average'),('branches',tables['branches'],6,0,'delta_average')]:
        for value in sorted({r[col_idx] for r in rows if r[col_idx] is not None}):
            for tolerance in (0,0.001,0.5):
                expected={r[id_idx] for r in rows if r[col_idx] is not None and value-tolerance<=r[col_idx]<=value+tolerance}
                column='rule_id' if table=='rules' else 'branch_id'
                actual={r[0] for r in con.execute(f'SELECT {column} FROM {table} WHERE {mass_column} BETWEEN ? AND ?',(value-tolerance,value+tolerance))}
                assert expected==actual,(table,value,tolerance,expected^actual)
                checks+=1
    con.commit();con.close();temp.replace(dest)
    manifest={'schema_version':'codex-index-1','source_logical_sha256':source_digest,'logical_sha256':logical,
              'database_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
              'counts':{k:len(v) for k,v in tables.items()},'direct_scan_comparisons':checks,
              'foreign_keys':'OK','dynamic_mass_policy':'instantiate_after_structural_applicability',
              'acceptance_scope':'index/scan equality; state-search acceptance is reported separately in validation/acceptance.json'}
    (ROOT/'index/manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    return manifest

if __name__=='__main__':print(json.dumps(build(),ensure_ascii=False))
