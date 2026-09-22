#!/usr/bin/env python3
"""Synthetic-only integration checks; explicit binaries, no shipped build/run."""
import argparse, ctypes, hashlib, json, os, pathlib, random, signal, subprocess, time

p = argparse.ArgumentParser()
p.add_argument('--binary', type=pathlib.Path, required=True)
p.add_argument('--baseline', type=pathlib.Path, required=True)
p.add_argument('--helper', type=pathlib.Path, required=True)
p.add_argument('--out', type=pathlib.Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
base_env = {k:v for k,v in os.environ.items() if not k.startswith('TESSERACT_')}
checks = 0

def check(ok, message):
    global checks
    assert ok, message
    checks += 1

def run(cmd, name, env=None, expected=0):
    e = dict(base_env); e.update(env or {})
    r = subprocess.run(list(map(str,cmd)), env=e, capture_output=True, timeout=40)
    (a.out/(name+'.log')).write_bytes(r.stdout+r.stderr)
    check(r.returncode == expected, (name,r.returncode,r.stderr.decode()[-2000:]))
    return r

def manifest(name, text):
    f = a.out/(name+'.tsv'); f.write_text('arm\tflags\n'+text); return f.resolve()

random.seed(71821)
genome = ''.join(random.choice('ACGT') for _ in range(12000))
rc = lambda s: s.translate(str.maketrans('ACGT','TGCA'))[::-1]
r1, r2 = a.out/'r1.fq', a.out/'r2.fq'
with r1.open('w') as x, r2.open('w') as y:
    for i in range(3000):
        pos = random.randrange(len(genome)-380)
        left, right = genome[pos:pos+151], rc(genome[pos+350-151:pos+350])
        x.write(f'@r{i}/1\n{left}\n+\n'+151*'I'+'\n')
        y.write(f'@r{i}/2\n{right}\n+\n'+151*'I'+'\n')
args = ['-1',r1.resolve(),'-2',r2.resolve(),'-k','31,55','-t','2','--no-html','--unitigs']

def asm(binary,out,extra=()): return [binary,*args,'-o',out,*extra]

def semantic(value):
    # Exact narrow resource/provenance exceptions; unknown fields still compare.
    excluded = {'started_at','command','total_seconds','peak_memory_bytes','seconds',
                'count_seconds','graph_seconds','simplify_seconds'}
    if isinstance(value,dict): return {k:semantic(v) for k,v in value.items() if k not in excluded}
    if isinstance(value,list): return [semantic(x) for x in value]
    return value

def same_outputs(x,y):
    for f in ['contigs.fasta','scaffolds.fasta','scaffolds.agp','assembly_graph.gfa','unitigs.fasta']:
        check((x/f).exists()==(y/f).exists(),f+' presence')
        if (x/f).exists(): check((x/f).read_bytes()==(y/f).read_bytes(),f+' bytes')
    check(semantic(json.loads((x/'report.json').read_text()))==semantic(json.loads((y/'report.json').read_text())), 'report semantic identity')

old,new = a.out/'old', a.out/'disabled'
run(asm(a.baseline,old),'old')
run(asm(a.binary,new),'disabled')
same_outputs(old,new)
check(not (new/'fork_child.json').exists(), 'disabled creates no fork artifact')

flags = ['TESSERACT_ROUTE_DISTANCE','TESSERACT_WEIGHTED_RESOLVER_COVERAGE',
         'TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE','TESSERACT_EXACT_READ_THREADS',
         'TESSERACT_OWNED_ANCHOR_PREFIX','TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT',
         'TESSERACT_OBSERVED_GRAPH_COVERAGE']
rows = 'first\t-\nowned\tTESSERACT_OWNED_ANCHOR_PREFIX=1\nobserved\tTESSERACT_OBSERVED_GRAPH_COVERAGE=1\nexplicit_off\t'+','.join(f+'=0' for f in flags)+'\nlast\t-\n'
m = manifest('batch',rows); out = a.out/'batch'
r = run(asm(a.binary,out),'batch',{'TESSERACT_DEV_FORK_BATCH':str(m),'TESSERACT_ROUTE_DISTANCE':'1','TESSERACT_OBSERVED_GRAPH_COVERAGE':'1'})
check(b'fork-batch complete' in r.stderr and b'  contigs      ' not in r.stderr,'parent batch completion has no assembly summary')
check(not (out/'contigs.fasta').exists() and not (out/'report.json').exists(), 'no parent assembly')
batch = json.loads((out/'fork_batch.json').read_text())
check(batch['status']=='complete' and len(batch['arms'])==5,'all arms completed')
check(batch['binary_sha256']==hashlib.sha256(a.binary.read_bytes()).hexdigest(),'executable SHA-256')
check(batch['manifest_sha256']==hashlib.sha256(m.read_bytes()).hexdigest(),'parsed manifest SHA-256')
check(batch['parent_tesseract_environment']['TESSERACT_ROUTE_DISTANCE']=='1','parent environment recorded')
check(batch['parent_tesseract_environment']['TESSERACT_OBSERVED_GRAPH_COVERAGE']=='1','inherited observed-depth environment recorded')
for arm in ['first','explicit_off','last']:
    same_outputs(new,out/arm)
for arm in ['first','owned','observed','explicit_off','last']:
    child=json.loads((out/arm/'fork_child.json').read_text())
    check(child['status']=='complete' and child['mode']=='fork-child','child completion provenance')
    check(child['actual_parent_invocation']==batch['actual_parent_invocation'],'actual parent command retained')
    check(child['effective_child_options']['outDir']==str((out/arm).resolve()),'isolated effective outdir')
    check('postgraph_child_only' in child['report_total_seconds_scope'],'child timer scope explicit')
    check('child_process_only' in child['report_peak_memory_bytes_scope'],'child peak-memory scope explicit')
    check((out/arm/'diagnostic.log').exists(),'isolated log')
    expected_observed = '1' if arm=='observed' else ('0' if arm=='explicit_off' else None)
    check(child['postgraph_flags']['TESSERACT_OBSERVED_GRAPH_COVERAGE']==expected_observed,'observed-depth effective child flag')
    check((b'[observed graph coverage]' in (out/arm/'diagnostic.log').read_bytes())==(arm=='observed'),'observed measurement runs only in enabled child')
check(batch['arms'][0]['pid']!=batch['arms'][1]['pid'],'distinct child PIDs')
check(all(x['sampled_parent_child_pss_sum_peak_bytes']>0 for x in batch['arms']), 'memory samples present')
owned = a.out/'owned_standalone'
run(asm(a.binary,owned),'owned_standalone',{'TESSERACT_OWNED_ANCHOR_PREFIX':'1'})
same_outputs(owned,out/'owned')
observed = a.out/'observed_standalone'
run(asm(a.binary,observed),'observed_standalone',{'TESSERACT_OBSERVED_GRAPH_COVERAGE':'1'})
same_outputs(observed,out/'observed')
# Depth is a real arm effect, while graph sequence/link construction stays fixed.
def graph_structure_and_depth(directory):
    structure=[];depth=[]
    for line in (directory/'assembly_graph.gfa').read_text().splitlines():
        fields=line.split('\t')
        if fields[0]=='S':
            structure.append(tuple(x for x in fields if not x.startswith('dp:f:')))
            depth.append(next(x for x in fields if x.startswith('dp:f:')))
        elif fields[0]=='L': structure.append(tuple(fields))
    return structure,depth
baseline_structure,baseline_depth=graph_structure_and_depth(new)
observed_structure,observed_depth=graph_structure_and_depth(out/'observed')
check(baseline_structure==observed_structure,'observed fork preserves final graph construction')
check(baseline_depth!=observed_depth,'observed fork actually changes graph depth')

# All effective option members are represented; future additions require review.
expected = set('libraries outDir kValues threads forcedCutoff trustCutoff minMaskRun minContigLen verbose correctReads resolveRepeats scaffold gapFill polish emitGfa emitHtml emitUnitigs mode minLinkSupport linkSupportPerX tieRatio bubbleCoverageLimit simplifyRounds polishPasses maxMemoryBytes qtrim ladderUnion ladderUnionMaxPresent ladderUnionMinLen dedupContained trimTerminalOverlap commandLine organism organismModelPath isPanelPath isSitesPath qcPath layout mapPolisher mapperDir userSetK userSetMinLink userSetTie userSetLinkPerX userSetBubble userSetRounds userSetPolishPasses'.split())
check(set(batch['effective_parent_options'])==expected,'complete AssemblyOptions inventory')
for name,text in [('upstream','bad\tTESSERACT_GAPCLOSE=1\n'),('ec','bad\tTESSERACT_EC_REQUIRE_UNIQUE_BEST=1\n'),('duplicate','a\t-\na\t-\n'),('path','../a\t-\n'),('value','a\tTESSERACT_ROUTE_DISTANCE=2\n')]:
    mf=manifest(name,text)
    r=run(asm(a.binary,a.out/name),'reject_'+name,{'TESSERACT_DEV_FORK_BATCH':str(mf)},expected=1)
    check(b'[1/7]' not in r.stderr,'invalid manifest rejected before reads')
for name,extra,env in [('no_resolve',['--no-resolve'],{}),('mapper',['--map-polish','bowtie2'],{}),('dump',[],{'TESSERACT_JOIN_DUMP':'/tmp/forbidden'}),('cluster',[],{'TESSERACT_PLASMID_CLUSTERS':'/tmp/forbidden'})]:
    r=run(asm(a.binary,a.out/name,extra),'reject_'+name,{'TESSERACT_DEV_FORK_BATCH':str(m),**env},expected=1)
    check(b'[1/7]' not in r.stderr,'incompatible mode rejected before reads')

# Hash independent known vectors; fixture uses the actual helper implementation.
for i, data in enumerate([b'',b'abc',bytes(range(256))*3+b'boundary',b'x'*65537]):
    f=a.out/f'hash{i}.bin';f.write_bytes(data)
    r=run([a.helper,'hash',f],'hash'+str(i))
    check(r.stdout.decode().strip()==hashlib.sha256(data).hexdigest(),'SHA-256 known/multiblock vector')

hm=manifest('helper','baseline\t-\nflagged\tTESSERACT_ROUTE_DISTANCE=1\nobserved\tTESSERACT_OBSERVED_GRAPH_COVERAGE=1\nlast\t-\n')
for action,expected_code in [('success',0),('exit',1),('signal',1),('thread',1),('sigchld_ignored',1),('sigchld_no_wait',1),('blocked',1)]:
    dest=a.out/('helper_'+action);dest.mkdir()
    run([a.helper,action,dest,r1,r2],action,{'TESSERACT_DEV_FORK_BATCH':str(hm),'TESSERACT_ROUTE_DISTANCE':'1','TESSERACT_OBSERVED_GRAPH_COVERAGE':'1'},expected_code)
    b=json.loads((dest/'fork_batch.json').read_text())
    check(b['status']==('complete' if action=='success' else 'failed'),'helper actual parent status')
    check(b['effective_parent_options']['maxMemoryBytes']==9007199254740993,'integer option provenance lossless')
    if action=='success':
        check(json.loads((dest/'baseline/report.json').read_text())['TESSERACT_ROUTE_DISTANCE'] is None,'actual child getenv cleared inherited flag')
        check(json.loads((dest/'flagged/report.json').read_text())['TESSERACT_ROUTE_DISTANCE']=='1','actual child getenv enabled')
        check(json.loads((dest/'last/report.json').read_text())['TESSERACT_ROUTE_DISTANCE'] is None,'later child independently cleared')
        check(json.loads((dest/'baseline/report.json').read_text())['TESSERACT_OBSERVED_GRAPH_COVERAGE'] is None,'actual child clears inherited observed flag')
        check(json.loads((dest/'observed/report.json').read_text())['TESSERACT_OBSERVED_GRAPH_COVERAGE']=='1','actual observed child getenv enabled')
        check(json.loads((dest/'last/report.json').read_text())['TESSERACT_OBSERVED_GRAPH_COVERAGE'] is None,'later child independently clears observed flag')
    if action=='exit': check(all(x['exit_code']==7 for x in b['arms']),'actual nonzero child exits propagate')
    if action=='signal': check(all(x['signal']==signal.SIGTERM for x in b['arms']),'actual child signals propagate')
    if action in ('thread','sigchld_ignored','sigchld_no_wait','blocked'):
        check(all(x['pid']==0 for x in b['arms']),'unsafe parent state never forks')

# Parent cancellation, inherited ignored disposition, and abrupt parent death.
# Become a temporary subreaper so the orphan child can be explicitly waited,
# rather than leaving a zombie to the container's unrelated init process.
libc=ctypes.CDLL(None, use_errno=True)
check(libc.prctl(36, 1, 0, 0, 0)==0,'test becomes child subreaper')
for action,death_signal in [('cancel',signal.SIGTERM),('ignored_cancel',signal.SIGTERM),('cancel',signal.SIGKILL)]:
    tag=action+'_'+str(death_signal)
    cancel=a.out/('helper_'+tag);cancel.mkdir()
    env=dict(base_env,TESSERACT_DEV_FORK_BATCH=str(hm))
    with (a.out/(tag+'.log')).open('wb') as log:
        parent=subprocess.Popen([str(a.helper),action,str(cancel),str(r1),str(r2)],env=env,stdout=log,stderr=log)
        deadline=time.monotonic()+10;pid=None
        while time.monotonic()<deadline:
            fp=cancel/'baseline/fork_child.json'
            if fp.exists():
                pid=json.loads(fp.read_text())['pid'];break
            time.sleep(.005)
        check(pid is not None,'cancellation child starts')
        parent.send_signal(death_signal)
        check(parent.wait(timeout=10)!=0,'parent cancellation fails batch')
        if death_signal==signal.SIGKILL:
            deadline=time.monotonic()+10;reaped=0;status=0
            while time.monotonic()<deadline:
                reaped,status=os.waitpid(pid,os.WNOHANG)
                if reaped: break
                time.sleep(.005)
            if not reaped:
                os.kill(pid,signal.SIGKILL);os.waitpid(pid,0)
            check(reaped==pid and os.WIFSIGNALED(status) and os.WTERMSIG(status)==signal.SIGTERM,'parent-death signal terminates orphan')
        check(not pathlib.Path(f'/proc/{pid}').exists(),'cancelled child reaped')
        b=json.loads((cancel/'fork_batch.json').read_text())
        if death_signal!=signal.SIGKILL:
            check(b['status']=='failed' and all(x['pid']==0 for x in b['arms'][1:]),'cancellation stops later forks')
        else:
            check(b['status']=='running','uncatchable death leaves honest incomplete metadata')
check(libc.prctl(36, 0, 0, 0, 0)==0,'test restores subreaper state')

# Claim must prevent a second batch before the first writes fork_batch.json,
# including manifests with disjoint arm IDs.
claim=a.out/'concurrent_claim';claim.mkdir()
other_manifest=manifest('other_claim','different_arm\t-\n')
with (a.out/'claim_holder.log').open('wb') as log:
    holder=subprocess.Popen([str(a.helper),'prepare_hold',str(claim),str(r1),str(r2)],env=dict(base_env,TESSERACT_DEV_FORK_BATCH=str(hm)),stdout=log,stderr=log)
    try:
        deadline=time.monotonic()+10
        while time.monotonic()<deadline and not (claim/'claim_ready').exists(): time.sleep(.005)
        check((claim/'claim_ready').exists() and not (claim/'fork_batch.json').exists(),'first batch exclusively claims root before preprocessing')
        r=run([a.helper,'success',claim,r1,r2],'concurrent_claim_reject',{'TESSERACT_DEV_FORK_BATCH':str(other_manifest)},expected=3)
        check(b'already claimed' in r.stderr,'concurrent root rejected by atomic claim')
    finally:
        holder.terminate();holder.wait(timeout=10)
check((claim/'.fork_batch_claim').exists(),'claim persists after interrupted preprocessing')
(a.out/'test_summary.json').write_text(json.dumps({'checks_passed':checks,'batch_memory_samples':[x['sampled_parent_child_pss_sum_peak_bytes'] for x in batch['arms']]},indent=2)+'\n')
print(f'{checks} checks passed; synthetic-only; no real libraries')
