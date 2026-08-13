# Apptainer / Singularity Container Architecture Guide

**Date**: 2026-08-12  
**Author**: Antigravity Assistant  
**Target Platform**: SGE HPC Cluster / Linux (x86_64)  
**Apptainer Version**: 1.3.2+  

---

## 1. Why Apptainer (.sif) is Superior to Conda on HPC

1. **Eliminates NFS I/O Metadata Thrashing**:
   - Conda environments store **50,000+ loose files** (shared libraries `.so`, Python `.pyc`, R headers) across network storage (`/data_lab_PGP/`).
   - When 20–50 parallel cluster jobs execute simultaneously, looking up loose files over NFS causes network latency and filesystem bottlenecks.
   - An Apptainer `.sif` image is a **single, immutable SquashFS file** read as contiguous blocks or cached directly into compute node RAM.

2. **100% Reproducible & Self-Contained**:
   - Bundles Ubuntu 22.04 base OS, glibc, Python 3.10, PyTorch, TensorFlow 2.15.1, SpliceAI 1.3.1, Pangolin, SPiP R runtime, `htslib`, `bcftools`, `samtools`, and `pysam`.
   - Immune to host OS package updates, library conflicts, or broken conda prefixes.

3. **Seamless Shared Filesystem Access**:
   - Mounts `/data_lab_PGP` and `/references` transparently via `-B /data_lab_PGP:/data_lab_PGP -B /references:/references`.

---

## 2. Container Specifications

- **Definition File**: [`resources/containers/annotation_pipeline.def`](../resources/containers/annotation_pipeline.def)
- **Target Image**: `/data_lab_PGP/resources/sif_images/annotation_pipeline.sif`
- **Build Script**: [`src/hpc/build_container.sh`](../src/hpc/build_container.sh)

---

## 3. How to Build the Container

### Option A: Build via SGE Cluster Job (Recommended)
On the cluster login node, submit the build job:

```bash
qsub src/hpc/build_container.sh
```

### Option B: Build Directly on Compute / Login Node
```bash
apptainer build --fakeroot /data_lab_PGP/resources/sif_images/annotation_pipeline.sif resources/containers/annotation_pipeline.def
```

---

## 4. Universal Container Execution Pattern

Once built, any pipeline component can be executed with a single line:

```bash
SIF="/data_lab_PGP/resources/sif_images/annotation_pipeline.sif"

# Run Python scripts (SpliceAI, Pangolin, Branchpointer, filter_variants):
apptainer exec -B /data_lab_PGP:/data_lab_PGP -B /references:/references "$SIF" python3 src/python/annotate_spliceai.py ...

# Run R scripts (SPiP, statsJPO, plotting):
apptainer exec -B /data_lab_PGP:/data_lab_PGP -B /references:/references "$SIF" Rscript src/R/plot_annotation_results.R ...

# Run Genomics binaries:
apptainer exec -B /data_lab_PGP:/data_lab_PGP -B /references:/references "$SIF" bcftools ...
```
