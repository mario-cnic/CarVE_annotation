# Nextflow scratch work directory

Recorded: 2026-10-08T13:26:19+02:00.

Mario selected scratch for new launches only. The `standard` profile in `nextflow.config` now sets `workDir` to `/mnt/scratch/mruizp/annotation_pipeline_new/work`. The wrapper forwards `-work-dir` overrides unchanged. `local_dev` is unchanged; published results stay under the existing project output directory.

No existing work, results, caches or running jobs were moved. `/mnt/scratch` is not mounted in this local session, so no destination directory was created and cluster accessibility is unverified. Verify the destination is writable and shared across execution nodes on the next cluster fixture, including VEP container execution.

For an older run stored in the project’s `work/`, launch from its original launch directory with the original arguments plus `-work-dir /data_lab_PGP/pipelines/annotation_pipeline_new/work -resume`. Resume needs both task outputs in that work directory and the corresponding `.nextflow/cache`; changing this default does not migrate either.

Local validation: cached Nextflow 24.04.2 successfully resolved both profiles. `standard` selects the scratch work directory and SGE; `local_dev` selects the local executor without the scratch override. Both retain the existing published output root. `git diff --check` passed. No workflow or SGE command was executed. Cluster Nextflow 25.10.2 and actual scratch task execution remain unverified.
