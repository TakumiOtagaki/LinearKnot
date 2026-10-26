#!/usr/bin/env python3
"""Scaffold selection and CLI contract; uses a controlled LP wrapper, no download."""
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'scripts/scaffold/make_scaffold.py'
spec = importlib.util.spec_from_file_location('make_scaffold', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

seq = list('A' * 22)
mea = '((((....))))..((....))'
for i,j in module.pairs_of(mea).items(): seq[i-1],seq[j-1]='G','C'
sequence=''.join(seq)
bpp = {(1,12):.9,(2,11):.9,(3,10):.2,(4,9):.9,(15,22):.8,(16,21):.8}
assert module.scaffold(sequence, mea, bpp, .5, 1) == '((.(....).))..((....))'
assert module.scaffold(sequence, mea, bpp, .5, 2) == '((........))..((....))'
assert module.scaffold(sequence, mea, bpp, .5, 3) == '.' * 22
assert module.scaffold(sequence, None, bpp) == '((........))..((....))'
assert module.scaffold(sequence, '.' * 22, bpp) == '.' * 22
assert module.scaffold('GGAAAAAACC', None, {(1,10):.7,(2,9):.7}) == '((......))'
assert module.scaffold('GGAAAAAACC', None, {(1,10):.7,(2,9):.699999}) == '.' * 10
assert module.scaffold('AAAAAA', None, {}) == '.' * 6

for bad in ({(1,10):.9,(1,9):.8}, {(1,10):.9,(2,10):.8}, {(1,8):.9,(3,10):.8}):
    try: module.scaffold('GGGGCCCCCC', None, bad)
    except ValueError: pass
    else: raise AssertionError('accepted conflicting selected pairs')
try: module.scaffold('AAAAAAAAAA', None, {(1,10):.9}, .7, 1)
except ValueError: pass
else: raise AssertionError('accepted noncanonical pair')
for bad in ('1 10 nan', '1 10 1.1', '0 10 .8', '1 11 .8', '1 10 .8\n1 10 .8', '1 10', 'garbage'):
    try: module.read_bpp(io.StringIO(bad), 10)
    except ValueError: pass
    else: raise AssertionError('accepted malformed BPP: '+bad)
assert module.read_bpp(io.StringIO('# header\n>seq\n1 10 .7\n'), 10) == {(1,10):.7}
assert list(module.read_sequences(io.StringIO('>x\nacgu\nacgu\n>y\nGGG\n'))) == [('x','acguacgu'),('y','GGG')]

with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
    tmp=Path(tmp); wrapper=tmp/'linearpartition'
    wrapper.write_text('#!'+sys.executable+'\n'+r'''
import json, os, pathlib, sys
args=sys.argv[1:]
with open(os.environ['MOCK_LP_LOG'],'a') as f:f.write(json.dumps(args)+'\n')
s=sys.stdin.read().strip()
if '-M' in args:print('((((....))))..........' if len(s)==22 else '.'*len(s))
with open(args[args.index('-r')+1],'w') as f:
    if len(s)==22:f.write('1 12 .9\n2 11 .9\n3 10 .2\n4 9 .9\n15 22 .8\n16 21 .8\n')
''')
    wrapper.chmod(0o755)
    log=tmp/'calls.jsonl';meta=tmp/'metadata.json'
    env=dict(os.environ,MOCK_LP_LOG=str(log))
    def run(*args, text=None):
        return subprocess.run([sys.executable,str(SCRIPT),'--linearpartition',str(wrapper),*args],
                              input=text if text is not None else '>rna\n'+sequence+'\n',
                              capture_output=True,text=True,env=env,cwd=ROOT)
    r=run('--metadata',str(meta));assert r.returncode==0,r.stderr
    lines=r.stdout.splitlines();assert lines[0].startswith('# linearknot_scaffold_metadata=')
    assert lines[1:]==['>rna',sequence,'((........))..((....))']
    config=json.loads(lines[0].split('=',1)[1]);saved=json.loads(meta.read_text())
    assert config['mode']=='bpp_threshold' and config['threshold']==.7 and config['min_helix']==2
    assert config['beam']==200 and config['dangles']==2 and config['gamma'] is None
    assert config['linearpartition']['commit'] is None  # Do not mislabel the outer repo's revision.
    assert config['linearpartition']['wrapper_sha256']==hashlib.sha256(wrapper.read_bytes()).hexdigest()
    assert saved['records'][0]['scaffold_sha256']==hashlib.sha256(lines[-1].encode()).hexdigest()
    call=json.loads(log.read_text().splitlines()[-1])
    assert '-M' not in call and '--gamma' not in call
    assert call[call.index('-b')+1]=='200' and call[call.index('-c')+1]=='0.7'
    assert call[call.index('--dangles')+1]=='2'
    old=run('--from-mea','--no-metadata-header');assert old.returncode==0,old.stderr
    assert old.stdout.splitlines()==['>rna',sequence,'((.(....).))..........']
    call=json.loads(log.read_text().splitlines()[-1])
    assert '-M' in call and '-c' not in call and call[call.index('-b')+1]=='100'
    assert call[call.index('--gamma')+1]=='3.0'
    explicit=run('--from-mea','--gamma','3','--threshold','0.5','--min-helix','1','--no-metadata-header')
    assert explicit.stdout==old.stdout
    dots=run(text='tttttt\n');assert dots.returncode==0,dots.stderr
    assert dots.stdout.splitlines()[1:]==['>seq1','UUUUUU','......']
    many=run(text='>one\n'+sequence+'\n>two\n'+sequence+'\n')
    assert many.returncode==0 and many.stdout.count('# linearknot_scaffold_metadata=')==1
    assert many.stdout.count('((........))..((....))')==2
    for bad in (['--threshold','nan'],['--threshold','inf'],['--threshold','1.1'],['--min-helix','0'],
                ['--beam','0'],['--gamma','3'],['--from-mea','--gamma','nan']):
        result=run(*bad);assert result.returncode!=0,bad
    bad=run(text='>bad\nACGN\n');assert bad.returncode!=0
print('test_make_scaffold=ok default_BPP legacy_MEA CLI_metadata malformed_BPP empty_and_batch')
