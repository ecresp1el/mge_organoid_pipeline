# Overnight URD Run Check

Date submitted: 2026-06-12

Two separate checkpointed/full URD jobs were submitted so they can run overnight and be checked tomorrow morning.

## Jobs

| Dataset | Job ID | Run label | Status at last check | Node | Purpose |
|---|---:|---|---|---|---|
| DIV30 | 51691451 | `div30_first_urd_paper_radial_glia_30k_checkpoint_v2` | Running | `gl3307` | Checkpointed 30k rerun with new stage/resource logs before any all-cell DIV30 run |
| DIV90 | 51691452 | `div90_urd_jia_lineage_full_v4_glia_tips` | Running | `gl3112` | Full retained-cell DIV90 run using frozen v4 Jia-lineage/glia-tip logic |

## Logs

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51691451.log
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51691452.log
```

## Output Roots

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_30k_checkpoint_v2/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips/
```

## First Things To Check

```bash
squeue -j 51691451,51691452
sacct -j 51691451,51691452 --format=JobID,JobName%32,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
```

Check these files if the jobs are still running or if a run stops unexpectedly:

```text
tables/slurm_resource_heartbeat.tsv
tables/div30_first_urd_stage_timings.tsv
lineage_tree_jia_endpoint_tips_v1/tables/lineage_tree_stage_timings.tsv
lineage_tree_jia_endpoint_tips_v1/tables/slurm_resource_heartbeat.tsv
```

## Interpretation

DIV90 full is the intended production-scale v4 run.

DIV30 is intentionally not all-cells yet. It is a checkpointed 30k rerun because the previous 30k job consumed substantial resources but did not leave the expected final `div30_first_urd_object.rds`. If this checkpointed 30k run completes cleanly and writes the expected stage/resource logs, then DIV30 can be promoted to all-cells next.
