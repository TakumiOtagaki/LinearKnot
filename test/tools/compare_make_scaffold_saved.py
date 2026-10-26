#!/usr/bin/env python3
"""Optional integration check against the paper checkout's frozen biological scaffolds."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'scripts/scaffold/make_scaffold.py'
PIN = 'b450fb3e63189073b68d385589035f992080aa3a'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--paper',type=Path,default=ROOT.parent/'LinearCPartyPaper/experiments')
    p.add_argument('--linearpartition',type=Path,required=True)
    p.add_argument('--python',default=sys.executable)
    p.add_argument('--legacy-helper',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    # Preserve the venv executable's directory: resolving its symlink would
    # otherwise select the system interpreter without python-gflags.
    env=dict(os.environ,PATH=str(Path(args.python).absolute().parent)+os.pathsep+os.environ.get('PATH',''),
             OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LC_ALL='C')
    selected=[('bprna-pilot-20261008',name) for name in
              ('bpRNA_RFAM_15825','bpRNA_RFAM_25128','bpRNA_CRW_30304','bpRNA_RFAM_27342')]
    selected += [('rnacentral-initial-20261008',name) for name in ('URS000233CDDE','URS00019626CE')]
    rows=[]
    spec=importlib.util.spec_from_file_location('make_scaffold',SCRIPT)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    for dataset,name in selected:
        folder=args.paper/'data'/dataset
        manifest=json.loads((folder/'dataset.json').read_text())
        saved=next(r for r in manifest['records'] if r['id']==name)
        source=folder/saved['input_file'];data=source.read_text().strip().splitlines()
        assert digest(source)==saved['input_sha256']
        fasta='\n'.join(data[:2])+'\n'
        sidecar=args.output/(name+'.json')
        command=[args.python,str(SCRIPT),'--linearpartition',str(args.linearpartition.resolve()),'--metadata',str(sidecar.resolve())]
        start=time.monotonic()
        result=subprocess.run(command,input=fasta,text=True,capture_output=True,env=env,cwd=ROOT,timeout=180,check=True)
        (args.output/(name+'.stdout')).write_text(result.stdout)
        (args.output/(name+'.stderr')).write_text(result.stderr)
        payload=[line for line in result.stdout.splitlines() if not line.startswith('#')]
        assert payload==data,(name,'saved scaffold differs')
        metadata=json.loads(sidecar.read_text())
        assert metadata['linearpartition']['commit']==PIN
        assert metadata['linearpartition']['tracked_source_modified'] is False
        assert (metadata['mode'],metadata['beam'],metadata['threshold'],metadata['min_helix'],metadata['dangles'])==('bpp_threshold',200,.7,2,2)
        with (folder/'scaffold_logs'/name/'pairs.tsv').open() as f:bpp=helper.read_bpp(f,len(data[1]))
        candidates={(i,j) for (i,j),prob in bpp.items() if prob>=.7}
        isolated={p for p in candidates if (p[0]-1,p[1]+1) not in candidates and (p[0]+1,p[1]-1) not in candidates}
        row={'id':name,'length':len(data[1]),'dataset':dataset,'input_sha256':digest(source),
             'saved_scaffold_sha256':saved['scaffold_sha256'],'actual_scaffold_sha256':metadata['records'][0]['scaffold_sha256'],
             'isolated_threshold_pairs':len(isolated),'threshold_pairs':len(candidates),'command':command,
             'elapsed_seconds':time.monotonic()-start,'match':True}
        assert row['saved_scaffold_sha256']==row['actual_scaffold_sha256']
        if args.legacy_helper and len(rows)<2:
            common=['--linearpartition',str(args.linearpartition.resolve())]
            old=subprocess.run([args.python,str(args.legacy_helper.resolve()),*common],input=fasta,text=True,capture_output=True,env=env,cwd=ROOT,timeout=120,check=True)
            new=subprocess.run([args.python,str(SCRIPT),*common,'--from-mea','--no-metadata-header'],input=fasta,text=True,capture_output=True,env=env,cwd=ROOT,timeout=120,check=True)
            assert new.stdout==old.stdout,(name,'legacy MEA output differs')
            row['legacy_mode_byte_identical']=True
        rows.append(row)
        print(name+' scaffold=identical elapsed=%.2fs'%row['elapsed_seconds'],flush=True)
    assert any(r['isolated_threshold_pairs'] for r in rows)
    assert any(r['threshold_pairs']==0 for r in rows)
    report={'linearpartition_commit':PIN,'helper_sha256':digest(SCRIPT),'wrapper_sha256':digest(args.linearpartition),
            'python':args.python,'cases':rows,'status':'passed'}
    (args.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('saved_scaffold_comparison=ok inputs='+str(len(rows)))


if __name__=='__main__':
    main()
