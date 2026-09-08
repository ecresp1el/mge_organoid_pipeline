"""Immutable, individually restartable analysis checkpoints and artifact checks.

An attempt is reusable only after all saved artifacts reopen successfully,
their hashes are sealed, and COMPLETE.json is published. CURRENT.json points
to a complete attempt. Failed attempts are retained, never overwritten. This
module has no clustering logic and cannot submit jobs.
"""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import subprocess
import time
import uuid
import numpy as np
import pandas as pd
from scipy import sparse
from .provenance import sha256, write_json


def reopen_artifact(path):
    """Reopen every supported saved format, returning its structural summary.

    Binary Annoy models require the adjacent PCA dimension and cell table.
    Domain alignment checks are performed separately by the stage validator.
    Large tabular files are streamed in chunks rather than loaded together.
    """
    path=Path(path);name=path.name
    if name.endswith('.npy'):
        value=np.load(path,mmap_mode='r',allow_pickle=False)
        return dict(format='npy',shape=list(value.shape),dtype=str(value.dtype))
    if name.endswith('.npz'):
        with np.load(path,allow_pickle=False) as archive:
            summary={k:dict(shape=list(archive[k].shape),dtype=str(archive[k].dtype)) for k in archive.files}
            is_sparse='format' in archive.files and 'shape' in archive.files
        if is_sparse:
            matrix=sparse.load_npz(path);matrix.check_format(full_check=True)
            if not np.isfinite(matrix.data).all():raise ValueError('Nonfinite sparse data: '+str(path))
            return dict(format='sparse_npz',shape=list(matrix.shape),nnz=int(matrix.nnz),dtype=str(matrix.dtype))
        return dict(format='npz',arrays=summary)
    if name.endswith('.json'):
        value=json.loads(path.read_text())
        return dict(format='json',root_type=type(value).__name__)
    if name.endswith('.jsonl'):
        count=0
        with path.open() as handle:
            for line in handle:
                json.loads(line);count+=1
        return dict(format='jsonl',records=count)
    if name.endswith(('.tsv','.tsv.gz')):
        count=0;columns=[]
        try:
            for chunk in pd.read_csv(path,sep='\t',chunksize=10000):
                count+=len(chunk);columns=list(chunk.columns)
        except pd.errors.EmptyDataError:
            return dict(format='tsv',rows=0,columns=[],explicitly_empty=True)
        return dict(format='tsv',rows=count,columns=columns)
    if name=='annoy.index':
        from annoy import AnnoyIndex
        model=pd.read_csv(path.parent/'pca_components_retained.tsv',sep='\t',index_col=0)
        index=AnnoyIndex(model.shape[1],'euclidean');index.load(str(path))
        count=index.get_n_items();index.unload()
        return dict(format='annoy',items=count,dimensions=model.shape[1])
    if path.suffix in ['.txt','.md','.log','.out','.err','.py','.R','.sbatch']:
        with path.open(errors='strict') as handle:
            for _ in handle:pass
        return dict(format='text')
    raise ValueError('No reopen validator for saved artifact '+str(path))


def artifact_manifest(directory):
    """Reopen and hash all attempt artifacts after the writer has closed them."""
    records=[]
    for path in sorted(Path(directory).rglob('*')):
        if not path.is_file() or path.name in ['ARTIFACTS.json','COMPLETE.json','FAILURE.json']:
            continue
        structure=reopen_artifact(path)
        records.append(dict(path=str(path.relative_to(directory)),bytes=path.stat().st_size,
                            sha256=sha256(path),structure=structure))
    return records


def complete_attempt(stage_dir, expected_contract=None, deep=True):
    """Return the verified current attempt or raise for absent/corrupt output."""
    stage_dir=Path(stage_dir)
    current=json.loads((stage_dir/'CURRENT.json').read_text())
    attempt=stage_dir/current['attempt']
    complete=json.loads((attempt/'COMPLETE.json').read_text())
    if expected_contract is not None and complete['contract']!=expected_contract:
        raise ValueError('Completed stage belongs to a different frozen contract')
    manifest=attempt/'ARTIFACTS.json'
    if sha256(manifest)!=complete['manifest_sha256'] or complete['manifest_sha256']!=current['manifest_sha256']:
        raise ValueError('Completion seal mismatch')
    records=json.loads(manifest.read_text())
    for record in records:
        path=attempt/record['path']
        if not path.is_file() or path.stat().st_size!=record['bytes'] or sha256(path)!=record['sha256']:
            raise ValueError('Artifact missing or changed: '+str(path))
        if deep and reopen_artifact(path)!=record['structure']:
            raise ValueError('Artifact structure changed: '+str(path))
    return attempt,complete


def checkpoint(stage_dir, contract, compute, validate):
    """Reuse verified output or run a new exclusive attempt and publish it.

    ``compute(attempt)`` writes only inside that new attempt. ``validate`` checks
    scientific dimensions/IDs in reopened outputs and returns JSON evidence.
    A corrupt completed attempt is not silently reused or deleted. A stale
    execution lock is archived only after its recorded process/job is inactive.
    """
    stage_dir=Path(stage_dir);stage_dir.mkdir(parents=True,exist_ok=True)
    if (stage_dir/'CURRENT.json').exists():
        try:
            attempt,_=complete_attempt(stage_dir,contract,deep=True)
            print(json.dumps(dict(checkpoint=str(stage_dir),action='REUSE_VERIFIED',attempt=str(attempt))),flush=True)
            return attempt
        except (ValueError,FileNotFoundError):
            # Preserve the failed seal for inspection, then make a new attempt.
            (stage_dir/'CURRENT.json').rename(stage_dir/('INVALID_CURRENT_'+uuid.uuid4().hex+'.json'))
    lock=stage_dir/'RUNNING.json'
    if lock.exists():
        owner=json.loads(lock.read_text());job=owner.get('slurm_job_id')
        if job:
            active=subprocess.check_output(['squeue','-h','-j',str(job)],text=True).strip()
            if active:raise RuntimeError('Checkpoint already owned by active job '+str(job))
        else:
            try:os.kill(owner['pid'],0)
            except ProcessLookupError:pass
            else:raise RuntimeError('Checkpoint already owned by live process')
        lock.rename(stage_dir/('STALE_LOCK_'+uuid.uuid4().hex+'.json'))
    descriptor=os.open(str(lock),os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    with os.fdopen(descriptor,'w') as handle:
        json.dump(dict(pid=os.getpid(),slurm_job_id=os.environ.get('SLURM_JOB_ID')),handle)
    name=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:8]
    attempt=stage_dir/'attempts'/name;attempt.mkdir(parents=True)
    start=time.monotonic();usage0=resource.getrusage(resource.RUSAGE_SELF)
    try:
        compute(attempt)
        evidence=validate(attempt)
        write_json(attempt/'VALIDATION.json',evidence)
        usage=resource.getrusage(resource.RUSAGE_SELF)
        write_json(attempt/'resources.json',dict(wall_seconds=time.monotonic()-start,
            cpu_seconds=usage.ru_utime+usage.ru_stime-usage0.ru_utime-usage0.ru_stime,
            peak_rss_gib=usage.ru_maxrss/1024**2,slurm_job_id=os.environ.get('SLURM_JOB_ID')))
        records=artifact_manifest(attempt)
        write_json(attempt/'ARTIFACTS.json',records)
        seal=sha256(attempt/'ARTIFACTS.json')
        write_json(attempt/'COMPLETE.json',dict(contract=contract,manifest_sha256=seal,
            artifacts=len(records),status='COMPLETE_VALIDATED',utc=datetime.now(timezone.utc).isoformat()))
        write_json(stage_dir/'CURRENT.json',dict(attempt=str(attempt.relative_to(stage_dir)),manifest_sha256=seal))
        # Independent read after publication: completion is never based on a
        # writer's in-memory objects alone.
        complete_attempt(stage_dir,contract,deep=True)
        return attempt
    except Exception as error:
        write_json(attempt/'FAILURE.json',dict(type=type(error).__name__,message=str(error)))
        raise
    finally:
        lock.unlink(missing_ok=True)
