#!/usr/bin/env python3
"""Integration regression: EC budget fails before trusting partial counts.

Usage: python3 tests/test_ec_memory_budget.py PATCHED_BINARY BASELINE_BINARY OUT
Both binaries remain untouched; all artifacts are written under OUT.
"""
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys

patched,baseline=map(lambda s:Path(s).resolve(),sys.argv[1:3])
out=Path(sys.argv[3]).resolve()
out.mkdir(parents=True,exist_ok=True)
rng=random.Random(49261)
genome=''.join(rng.choice('ACGT') for _ in range(12000))
rev=str.maketrans('ACGT','TGCA')
r1,r2=out/'reads1.fq',out/'reads2.fq'
with r1.open('w') as a,r2.open('w') as b:
 for i in range(1200):
  pos=rng.randrange(len(genome)-350+1)
  x=genome[pos:pos+150];y=genome[pos+200:pos+350].translate(rev)[::-1]
  a.write(f'@pair{i}/1\n{x}\n+\n'+150*'I'+'\n')
  b.write(f'@pair{i}/2\n{y}\n+\n'+150*'I'+'\n')
env={k:v for k,v in os.environ.items() if not k.startswith('TESSERACT_')}
results={}
for label,binary,memory in [('budget_low',patched,'.001'),('enough',patched,'1'),('more',patched,'2'),('baseline',baseline,'1')]:
 target=out/label
 if target.exists():raise RuntimeError(f'refusing to reuse existing test output {target}')
 cmd=[str(binary),'-1',str(r1),'-2',str(r2),'-o',str(target),'-k','21,33,55','-t','1','--no-html','--no-gfa','--max-memory',memory]
 run=subprocess.run(cmd,text=True,capture_output=True,env=env,timeout=120)
 (out/f'{label}.log').write_text(run.stdout+run.stderr)
 if label=='budget_low':
  assert run.returncode!=0,'low memory must fail'
  assert 'while counting 21-mers for read error correction' in run.stderr,run.stderr
  assert not (target/'report.json').exists(),'failed run cannot publish completed report'
  assert not (target/'contigs.fasta').exists(),'failed run cannot publish partial contigs'
  results[label]={'exit_code':run.returncode,'error_stage':'EC counting','completed_outputs':False}
 else:
  assert run.returncode==0,run.stderr
  report=json.loads((target/'report.json').read_text())
  fasta=(target/'contigs.fasta').read_bytes()
  assert fasta.startswith(b'>') and report['error_correction']['run'],'successful EC assembly must be nonempty'
  results[label]={'exit_code':0,'contigs_sha256':hashlib.sha256(fasta).hexdigest()}
assert results['enough']['contigs_sha256']==results['more']['contigs_sha256']==results['baseline']['contigs_sha256'],'available memory must not change assembled sequence'
(out/'result.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
