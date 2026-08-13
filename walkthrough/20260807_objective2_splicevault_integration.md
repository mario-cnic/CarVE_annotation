# Objective 2 Completion Walkthrough: SpliceVault Empirical RNA Integration

**Date**: 2026-08-07  
**Author**: Bioinformatics Pipeline Team  
**Status**: ✅ COMPLETED  

---

## 1. Summary of Accomplishments

All tasks for **Objective 2: SpliceVault Empirical RNA Integration** have been completed and verified across the pipeline:

### Key Enhancements & Validations:

1. **Dataset & Plugin Verification** ([annotate_vep_vars.sh](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_vep_vars.sh)):
   - Verified Ensembl VEP 111 SpliceVault plugin configuration:
     `--plugin SpliceVault,file=/data_lab_PGP/resources/annotation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz`
   - Verified dataset file integrity: `SpliceVault_data_GRCh38.tsv.gz` (883,907,365 bytes) and index `.tbi` (747,762 bytes).

2. **Column Registry Completeness** ([all_but_old_gnomad_vep_cols.txt](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/all_but_old_gnomad_vep_cols.txt)):
   - Confirmed registration of all 7 VEP SpliceVault output columns:
     - `SpliceVault_top_events`
     - `SpliceVault_site_sample_count`
     - `SpliceVault_site_max_depth`
     - `SpliceVault_out_of_frame_events`
     - `SpliceVault_site_pos`
     - `SpliceVault_site_type`
     - `SpliceVault_SpliceAI_delta`

3. **Status Contract & Decoded Event Extraction** ([filter_variants.py](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py)):
   - Enhanced `parse_splicevault()` in `filter_variants.py` to create `SpliceVault_status`:
     - `aberrant_event_detected`: Witnessed aberrant splicing events (exon skipping, cryptic donor/acceptor) in $>300,000$ RNA-seq sample cohort.
     - `no_events_found`: Position/splice site covered in dataset but zero aberrant splicing witnessed.
     - `not_covered`: Variant outside dataset reach.
   - Decodes rank, event type, relative genomic offset, sample percentage support, and reading frame changes into human-readable descriptions (`SpliceVault_Predictions_Decoded`).

---

## 2. Verification & Audit Results

Audited pipeline execution across 380,292 variants across 7 cardiac genes (`BAG3`, `DSP`, `FLNC`, `LMNA`, `MYBPC3`, `PKP2`, `TTN`):
- **SpliceVault Column Schema Consistency**: 100%
- **Status Contract Validity**: All variants cleanly classified as `aberrant_event_detected`, `no_events_found`, or `not_covered`.
