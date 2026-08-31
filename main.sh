#!/bin/bash
# Universal Genomic Variant Annotation, Filtering, and Visualization Pipeline
# Master Orchestration Script (SGE Cluster Execution)
# Features: State-Aware Checkpointing, Modular Predictor Skipping, Dynamic SGE Dependency Chaining, Built-in Run Quality Auditor

set -euo pipefail

# Default Parameters
RAW_MASTER_DIR=""
GENE_SELECTED=""
INPUT_BUILD="GRCh38"
OUTPUT_FORMAT="pq"
COLUMN_PARAM=""
SKIP_GENES=""
SLEEP_TIME=5
RUN_NAME=""
AUDIT_RUN=""
audit_only=false
is_test_mode=false

# Customizable Filtering Parameters
MAX_AF=""
MIN_REVEL=""
MIN_ALPHAMISSENSE=""
MIN_SPIP=""
MIN_CADD=""
CONSEQUENCES=""

# Pipeline Checkpoint & Overwrite Flags (Default: Smart Resume / Skip Completed)
overwrite_all=false
force_varconv=false;    skip_varconv=false
force_spip=false;       skip_spip=false
force_vep=false;        skip_vep=false
force_pangolin=false;   skip_pangolin=false
force_spliceai=false;   skip_spliceai=false
force_branchpoint=false; skip_branchpoint=false
force_merge=false;      skip_merge=false
force_vcf2parsed=false; skip_vcf2parsed=false
force_reports=false;    skip_reports=false

# Parse Command Line Arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --raw-dir)
            RAW_MASTER_DIR="$2"
            shift 2
            ;;
        --run-name)
            RUN_NAME="$2"
            shift 2
            ;;
        --gene)
            GENE_SELECTED="$2"
            shift 2
            ;;
        --build)
            INPUT_BUILD="$2"
            shift 2
            ;;
        --output-format)
            OUTPUT_FORMAT="${2#.}"
            shift 2
            ;;
        --column-param)
            COLUMN_PARAM="$2"
            shift 2
            ;;
        --max-af)
            MAX_AF="$2"
            shift 2
            ;;
        --min-revel)
            MIN_REVEL="$2"
            shift 2
            ;;
        --min-alphamissense)
            MIN_ALPHAMISSENSE="$2"
            shift 2
            ;;
        --min-spip)
            MIN_SPIP="$2"
            shift 2
            ;;
        --min-cadd)
            MIN_CADD="$2"
            shift 2
            ;;
        --consequences)
            CONSEQUENCES="$2"
            shift 2
            ;;
        --sleep-time)
            SLEEP_TIME="$2"
            shift 2
            ;;
        --skip-genes)
            SKIP_GENES="$2"
            shift 2
            ;;
        --audit-run)
            AUDIT_RUN="$2"
            shift 2
            ;;
        --audit-only)
            audit_only=true
            shift 1
            ;;
        --test)
            is_test_mode=true
            RAW_MASTER_DIR="test_data/raw_vcfs"
            shift 1
            ;;
        # Overall & Granular Force / Skip Flags
        --overwrite-all)
            overwrite_all=true
            shift 1
            ;;
        --force-varconv)
            force_varconv=true
            shift 1
            ;;
        --skip-varconv)
            skip_varconv=true
            shift 1
            ;;
        --force-spip)
            force_spip=true
            shift 1
            ;;
        --skip-spip|--skip-spip-annotation)
            skip_spip=true
            shift 1
            ;;
        --force-vep)
            force_vep=true
            shift 1
            ;;
        --skip-vep|--skip-vep-annotation)
            skip_vep=true
            shift 1
            ;;
        --force-pangolin)
            force_pangolin=true
            shift 1
            ;;
        --skip-pangolin)
            skip_pangolin=true
            shift 1
            ;;
        --force-spliceai)
            force_spliceai=true
            shift 1
            ;;
        --skip-spliceai)
            skip_spliceai=true
            shift 1
            ;;
        --force-branchpoint)
            force_branchpoint=true
            shift 1
            ;;
        --skip-branchpoint)
            skip_branchpoint=true
            shift 1
            ;;
        --force-merge)
            force_merge=true
            shift 1
            ;;
        --skip-merge)
            skip_merge=true
            shift 1
            ;;
        --force-vcf2parsed)
            force_vcf2parsed=true
            shift 1
            ;;
        --skip-vcf2parsed)
            skip_vcf2parsed=true
            shift 1
            ;;
        --force-reports)
            force_reports=true
            shift 1
            ;;
        --skip-reports)
            skip_reports=true
            shift 1
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

PYTHON_EXE="/home/mruizp/apps/miniforge3/envs/datasci/bin/python3"
if [ ! -x "$PYTHON_EXE" ]; then
    PYTHON_EXE="/home/mruizp/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python3"
fi
if [ ! -x "$PYTHON_EXE" ]; then
    PYTHON_EXE="/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python3"
fi
if [ ! -x "$PYTHON_EXE" ]; then
    PYTHON_EXE="python3"
fi

R_EXE="/data_lab_PGP/shared/utils/conda_envs/datasci/bin/Rscript"
if [ ! -x "$R_EXE" ]; then
    R_EXE="Rscript"
fi

# Built-in Direct Run Quality Auditor Invocation
if [ -n "$AUDIT_RUN" ] || [ "$audit_only" = true ]; then
    target_audit="${AUDIT_RUN:-}"
    if [ -z "$target_audit" ] && [ -n "$RUN_NAME" ]; then
        target_audit="RUNS/${RUN_NAME}"
    fi
    if [ -z "$target_audit" ] && [ -n "$RAW_MASTER_DIR" ]; then
        target_audit="RUNS/$(basename "$RAW_MASTER_DIR")"
    fi
    echo "============================================================================"
    echo " Executing Master Pipeline Quality & Completeness Auditor"
    echo "============================================================================"
    $PYTHON_EXE src/python/audit_run_results.py \
        ${target_audit:+--run-dir "$target_audit"} \
        ${RAW_MASTER_DIR:+--raw-dir "$RAW_MASTER_DIR"}
    exit 0
fi

if [ -z "$RAW_MASTER_DIR" ]; then
    echo "============================================================================"
    echo " Universal Genomic Variant Annotation & Visualization Pipeline "
    echo "============================================================================"
    echo "Usage: $0 --raw-dir <DIR> [options]"
    echo ""
    echo "Options:"
    echo "  --raw-dir <DIR>            Directory containing input variant files (.xlsx, .csv, .tsv, .pq, .vcf)"
    echo "  --run-name <NAME>          Custom run folder name inside RUNS/ (default: folder name of --raw-dir)"
    echo "  --gene <GENE>              Filter to specific gene(s)"
    echo "  --skip-genes <LIST>        Comma-separated list of genes to skip"
    echo "  --build <BUILD>            Input genomic assembly build: hg19 or GRCh38 (default: GRCh38)"
    echo "  --output-format <FMT>      Output table format: pq, tsv, xlsx (default: pq)"
    echo "  --max-af <FLOAT>           Filter max gnomAD allele frequency threshold"
    echo "  --min-revel <FLOAT>        Filter min REVEL score threshold"
    echo "  --min-alphamissense <FLOAT> Filter min AlphaMissense score threshold"
    echo "  --min-spip <FLOAT>         Filter min SPiP splice score threshold"
    echo "  --min-cadd <FLOAT>         Filter min CADD phred score threshold"
    echo "  --consequences <LIST>      Comma-separated list of target VEP consequences"
    echo ""
    echo "Quality Audit & Standalone Diagnostics:"
    echo "  --test                     Run end-to-end test suite on test dataset into test_data/test_run/"
    echo "  --audit-run <DIR/NAME>     Audit quality, completeness %, and schemas of a run directory"
    echo "  --audit-only               Run only the audit report on the specified run without submitting jobs"
    echo ""
    echo "Modular Execution & Checkpoint Overrides:"
    echo "  --overwrite-all            Re-run all intermediate steps and overwrite results"
    echo "  --force-spliceai           Force re-run SpliceAI (-D 10000) regardless of existing output"
    echo "  --force-vep                Force re-run VEP annotation"
    echo "  --force-pangolin           Force re-run Pangolin splice predictor"
    echo "  --force-spip               Force re-run SPiP predictor"
    echo "  --force-branchpoint        Force re-run Branchpoint predictor pair"
    echo "  --force-merge              Force re-run multi-predictor merge"
    echo "  --force-vcf2parsed         Force re-run VCF to table parsing & filtering"
    echo "  --force-reports            Force re-generate interactive and clinical reports"
    echo "============================================================================"
    exit 1
fi

if [ ! -d "$RAW_MASTER_DIR" ]; then
    echo "============================================================================"
    echo " Error: Input directory does not exist: '$RAW_MASTER_DIR'"
    echo " Please provide a valid directory path containing your input variant files."
    echo "============================================================================"
    exit 1
fi

COUNT_FILES=$(find "$RAW_MASTER_DIR" -maxdepth 1 -type f \( -name "*.xlsx" -o -name "*.csv" -o -name "*.tsv" -o -name "*.pq" -o -name "*.parquet" -o -name "*.vcf*" \) | wc -l)
if [ "$COUNT_FILES" -eq 0 ]; then
    echo "============================================================================"
    echo " Error: No supported variant files found in '$RAW_MASTER_DIR'"
    echo " Supported extensions: .xlsx, .csv, .tsv, .pq, .parquet, .vcf, .vcf.gz"
    echo "============================================================================"
    exit 1
fi

if [ -n "$GENE_SELECTED" ]; then
    GENE_SELECTED=$(echo "$GENE_SELECTED" | tr '[:lower:]' '[:upper:]')
fi

if [ -n "$SKIP_GENES" ]; then
    SKIP_GENES=$(echo "$SKIP_GENES" | tr '[:lower:]' '[:upper:]')
fi

if [ "$is_test_mode" = true ]; then
    job_pfx="t_"
    if [ -z "$RUN_NAME" ]; then
        RUN_NAME="test_run"
    fi
    RUN_BASE_DIR="test_data/${RUN_NAME}"
else
    job_pfx=""
    RAW_FOLDER=$(basename "$RAW_MASTER_DIR")
    if [ -z "$RUN_NAME" ]; then
        RUN_NAME="$RAW_FOLDER"
    fi
    RUN_BASE_DIR="RUNS/${RUN_NAME}"
fi
TMP_MASTER_DIR="${RUN_BASE_DIR}/_tmp"
ANNOTATION_MASTER_DIR="${RUN_BASE_DIR}/annotation"
RESULTS_MASTER_DIR="${RUN_BASE_DIR}/results"
FILTERED_MASTER_DIR="${RUN_BASE_DIR}/filtered"
PLOTS_MASTER_DIR="${RUN_BASE_DIR}/plots"
REPORTS_MASTER_DIR="${RUN_BASE_DIR}/reports"
ERROR_LOG_DIR="${RUN_BASE_DIR}/_log"

GENES_BED_FILE="resources/cardio_genes_loc.bed"
GENE_TRANSCRIPT_MAPPING="resources/gene_transcript_mapping.txt"

mkdir -p "$RUN_BASE_DIR" "$TMP_MASTER_DIR" "$ANNOTATION_MASTER_DIR" "$RESULTS_MASTER_DIR" "$FILTERED_MASTER_DIR" "$PLOTS_MASTER_DIR" "$REPORTS_MASTER_DIR" "$ERROR_LOG_DIR"

# Helper: check non-empty existing file
is_valid_file() {
    [ -f "$1" ] && [ -s "$1" ]
}

echo "Resolving transcript mappings for target genes..."
$PYTHON_EXE src/python/query_new_transcripts.py "$RAW_MASTER_DIR" --mapping-file "$GENE_TRANSCRIPT_MAPPING" --auto-append || true

all_terminal_jobs=()

for input_file in "$RAW_MASTER_DIR"/*; do
    [ -f "$input_file" ] || continue
    ext="${input_file##*.}"
    if [[ "$ext" != "xlsx" && "$ext" != "pq" && "$ext" != "parquet" && "$ext" != "csv" && "$ext" != "tsv" && "$ext" != "vcf" && "$ext" != "gz" ]]; then
        continue
    fi

    base_name=$(basename "$input_file" | cut -d'.' -f1)
    gene_name=$(echo "$base_name" | cut -d'_' -f1 | tr '[:lower:]' '[:upper:]')

    if [ -n "$GENE_SELECTED" ] && [[ ! ",$GENE_SELECTED," == *",$gene_name,"* ]]; then
        continue
    fi

    if [ -n "$SKIP_GENES" ] && [[ ",$SKIP_GENES," == *",$gene_name,"* ]]; then
        echo "===================================================="
        echo " Skipping gene: $gene_name (specified in --skip-genes)"
        echo "===================================================="
        continue
    fi

    echo "----------------------------------------------------"
    echo " Submitting Annotation Pipeline Jobs for: $gene_name"
    echo "----------------------------------------------------"

    mkdir -p "$ERROR_LOG_DIR/$gene_name"
    final_output_file="$RESULTS_MASTER_DIR/${gene_name}.parsed.clean.${OUTPUT_FORMAT}"

    Transcript=$(grep -w "$gene_name" "$GENE_TRANSCRIPT_MAPPING" | cut -d',' -f3 || true)
    if [ -z "$Transcript" ]; then
        Transcript="UNKNOWN"
    fi

    # Step 0: Ensure BED coordinates exist
    if ! grep -q "$gene_name" "$GENES_BED_FILE"; then
        bash src/hpc/gene_coords.sh "$GENE_TRANSCRIPT_MAPPING" "$GENES_BED_FILE" "$gene_name"
    fi

    # Step 1: Variant Converter (universal_variant_converter.py)
    vcf_file="$TMP_MASTER_DIR/${gene_name}.vcf"
    vcf_gz="${vcf_file}.gz"
    var_conv_job=""

    if [ "$skip_varconv" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_varconv" = true ] || ! is_valid_file "$vcf_gz"; then
            var_conv_job=$(qsub -N "${job_pfx}varconv_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 1 \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.varconv.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.varconv.err" \
                -b y bash src/hpc/variant_converter.sh "$input_file" "$vcf_file" --build "$INPUT_BUILD" "$COLUMN_PARAM" | awk '{print $3}')
            echo "  [1/8] Submitted variant_converter job: $var_conv_job"
        else
            echo "  [1/8] Skipping variant_converter (already completed: $vcf_gz)"
        fi
    fi

    # Base dependency for parallel predictors (if Step 1 was submitted, wait for it)
    pred_hold_flag=""
    if [ -n "$var_conv_job" ]; then
        pred_hold_flag="-hold_jid $var_conv_job"
    fi

    active_predictor_jobs=()

    # Step 2.1: SPiP Annotation (direct from VCF.gz)
    spip_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annSPiP.vcf"
    spip_vcf_gz="${spip_vcf}.gz"
    spip_job=""
    if [ "$skip_spip" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_spip" = true ] || (! is_valid_file "$spip_vcf_gz" && ! is_valid_file "$spip_vcf"); then
            spip_job=$(qsub -N "${job_pfx}spip_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 4 \
                $pred_hold_flag \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spip.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spip.err" \
                -b y bash src/hpc/annotate_spip_vars.sh "$vcf_gz" "$spip_vcf" | awk '{print $3}')
            echo "  [2/8] Submitted SPiP annotation job: $spip_job"
            active_predictor_jobs+=("$spip_job")
        else
            echo "  [2/8] Skipping SPiP annotation (already completed)"
        fi
    fi

    # Step 2.2: VEP Annotation (direct from VCF.gz)
    vep_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annVEP.vcf.gz"
    vep_job=""
    if [ "$skip_vep" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_vep" = true ] || ! is_valid_file "$vep_vcf"; then
            vep_job=$(qsub -N "${job_pfx}vep_${gene_name}" -P BIGN -A PGP -l h_vmem=15G -pe smp 4 \
                $pred_hold_flag \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vep.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vep.err" \
                -b y bash src/hpc/annotate_vep_vars.sh "$vcf_gz" "$vep_vcf" | awk '{print $3}')
            echo "  [3/8] Submitted VEP annotation job: $vep_job"
            active_predictor_jobs+=("$vep_job")
        else
            echo "  [3/8] Skipping VEP annotation (already completed)"
        fi
    fi

    # Step 2.3: Pangolin Splicing Predictor (direct from VCF.gz)
    pangolin_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annPangolin.vcf.gz"
    pangolin_job=""
    if [ "$skip_pangolin" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_pangolin" = true ] || ! is_valid_file "$pangolin_vcf"; then
            pangolin_job=$(qsub -N "${job_pfx}pangolin_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 4 \
                $pred_hold_flag \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.pangolin.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.pangolin.err" \
                -b y bash src/hpc/annotate_pangolin_vars.sh "$vcf_gz" "$pangolin_vcf" | awk '{print $3}')
            echo "  [4/8] Submitted Pangolin annotation job: $pangolin_job"
            active_predictor_jobs+=("$pangolin_job")
        else
            echo "  [4/8] Skipping Pangolin annotation (already completed)"
        fi
    fi

    # Step 2.4: SpliceAI Predictor (Local -D 10000 with Parallel VCF Chunking)
    spliceai_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annSpliceAI.vcf.gz"
    spliceai_job=""
    if [ "$skip_spliceai" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_spliceai" = true ] || ! is_valid_file "$spliceai_vcf"; then
            spliceai_job=$(qsub -N "${job_pfx}spliceai_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 4 \
                $pred_hold_flag \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spliceai.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spliceai.err" \
                -b y bash src/hpc/annotate_spliceai_vars.sh "$vcf_gz" "$spliceai_vcf" | awk '{print $3}')
            echo "  [5/8] Submitted SpliceAI (-D 10000) annotation job: $spliceai_job"
            active_predictor_jobs+=("$spliceai_job")
        else
            echo "  [5/8] Skipping SpliceAI annotation (already completed: $spliceai_vcf)"
        fi
    fi

    # Step 2.5: Branchpoint Predictor Pair (Branchpointer + LaBranchoR)
    branchpoint_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annBranchpoint.vcf.gz"
    branchpoint_job=""
    if [ "$skip_branchpoint" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_branchpoint" = true ] || ! is_valid_file "$branchpoint_vcf"; then
            branchpoint_job=$(qsub -N "${job_pfx}branchpoint_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 2 \
                $pred_hold_flag \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.branchpoint.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.branchpoint.err" \
                -b y bash src/hpc/annotate_branchpointer_vars.sh "$vcf_gz" "$branchpoint_vcf" | awk '{print $3}')
            echo "  [6/8] Submitted Branchpoint annotation job: $branchpoint_job"
            active_predictor_jobs+=("$branchpoint_job")
        else
            echo "  [6/8] Skipping Branchpoint annotation (already completed)"
        fi
    fi

    # Step 3: Merge all Predictor Annotations into Final Annotated VCF
    final_annotated_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annotated.vcf.gz"
    large_sv_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.large_svs.vcf.gz"
    merge_ann_job=""

    merge_hold_flag=""
    if [ ${#active_predictor_jobs[@]} -gt 0 ]; then
        joined_pred_jobs=$(IFS=,; echo "${active_predictor_jobs[*]}")
        merge_hold_flag="-hold_jid $joined_pred_jobs"
    fi

    if [ "$skip_merge" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_merge" = true ] || [ ${#active_predictor_jobs[@]} -gt 0 ] || ! is_valid_file "$final_annotated_vcf"; then
            merge_ann_job=$(qsub -N "${job_pfx}merge_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
                $merge_hold_flag \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.merge.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.merge.err" \
                -b y bash src/hpc/merge_vep_spip.sh "$spip_vcf.gz" "$vep_vcf" "$final_annotated_vcf" "$large_sv_vcf" "$pangolin_vcf" "$spliceai_vcf" "$branchpoint_vcf" | awk '{print $3}')
            echo "  [7/8] Submitted multi-predictor merge job: $merge_ann_job"
        else
            echo "  [7/8] Skipping multi-predictor merge (already completed)"
        fi
    fi

    # Step 4: Direct VCF to Final Clean Table (Parquet/Excel/TSV)
    vcf2parsed_job=""
    vcf2parsed_hold_flag=""
    if [ -n "$merge_ann_job" ]; then
        vcf2parsed_hold_flag="-hold_jid $merge_ann_job"
    fi

    if [ "$skip_vcf2parsed" = false ]; then
        if [ "$overwrite_all" = true ] || [ "$force_vcf2parsed" = true ] || [ -n "$merge_ann_job" ] || ! is_valid_file "$final_output_file"; then
            vcf2parsed_job=$(qsub -N "${job_pfx}vcf2parsed_${gene_name}" -P BIGN -A PGP -l h_vmem=80G -pe smp 1 \
                $vcf2parsed_hold_flag \
                -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vcf2parsed.out" \
                -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vcf2parsed.err" \
                -b y bash src/hpc/vcf2parsed.sh "$final_annotated_vcf" "$final_output_file" "$Transcript" | awk '{print $3}')
            echo "  [8/8] Submitted direct vcf2parsed job: $vcf2parsed_job"
        else
            echo "  [8/8] Skipping direct vcf2parsed (already completed: $final_output_file)"
        fi
    fi

    # Step 5: Downstream Filtering & Reports Generation
    downstream_hold_flag=""
    if [ -n "$vcf2parsed_job" ]; then
        downstream_hold_flag="-hold_jid $vcf2parsed_job"
    fi

    if [ "$skip_reports" = false ]; then
        filter_flags="--input $final_output_file --output-dir $FILTERED_MASTER_DIR"
        [ -n "$MAX_AF" ] && filter_flags="$filter_flags --max-af $MAX_AF"
        [ -n "$MIN_REVEL" ] && filter_flags="$filter_flags --min-revel $MIN_REVEL"
        [ -n "$MIN_ALPHAMISSENSE" ] && filter_flags="$filter_flags --min-alphamissense $MIN_ALPHAMISSENSE"
        [ -n "$MIN_SPIP" ] && filter_flags="$filter_flags --min-spip $MIN_SPIP"
        [ -n "$MIN_CADD" ] && filter_flags="$filter_flags --min-cadd $MIN_CADD"
        [ -n "$CONSEQUENCES" ] && filter_flags="$filter_flags --consequences \"$CONSEQUENCES\""

        filter_job=$(qsub -N "${job_pfx}filter_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
            $downstream_hold_flag \
            -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.filter.out" \
            -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.filter.err" \
            -b y bash src/hpc/run_python_hpc.sh "$PYTHON_EXE" src/python/filter_and_summarize.py $filter_flags | awk '{print $3}')

        plot_job=$(qsub -N "${job_pfx}plot_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
            $downstream_hold_flag \
            -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.plot.out" \
            -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.plot.err" \
            -b y $R_EXE src/R/plot_annotation_results.R --input "$final_output_file" --output-dir "$PLOTS_MASTER_DIR" | awk '{print $3}')

        report_job=$(qsub -N "${job_pfx}report_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
            $downstream_hold_flag \
            -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.report.out" \
            -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.report.err" \
            -b y bash src/hpc/run_python_hpc.sh "$PYTHON_EXE" src/python/generate_interactive_report.py --input "$final_output_file" --output "$REPORTS_MASTER_DIR/${gene_name}_interactive_dashboard.html" | awk '{print $3}')

        clinical_report_job=$(qsub -N "${job_pfx}clinrep_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
            $downstream_hold_flag \
            -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.clinrep.out" \
            -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.clinrep.err" \
            -b y bash src/hpc/run_python_hpc.sh "$PYTHON_EXE" src/python/generate_clinical_prioritization_report.py --input "$final_output_file" --output "$REPORTS_MASTER_DIR/${gene_name}_clinical_prioritization_report.html" | awk '{print $3}')

        [ -n "$clinical_report_job" ] && all_terminal_jobs+=("$clinical_report_job")
        [ -n "$report_job" ] && all_terminal_jobs+=("$report_job")
        [ -n "$plot_job" ] && all_terminal_jobs+=("$plot_job")
        [ -n "$filter_job" ] && all_terminal_jobs+=("$filter_job")
    elif [ -n "$vcf2parsed_job" ]; then
        all_terminal_jobs+=("$vcf2parsed_job")
    fi

    sleep $SLEEP_TIME
done

# Step 6: Run-Level Automated Quality & Completeness Auditor Job
if [ ${#all_terminal_jobs[@]} -gt 0 ]; then
    audit_hold_ids=$(IFS=,; echo "${all_terminal_jobs[*]}")
    audit_job=$(qsub -N "${job_pfx}audit_${RUN_NAME}" -P BIGN -A PGP -l h_vmem=15G -pe smp 1 \
        -hold_jid "$audit_hold_ids" \
        -o "$ERROR_LOG_DIR/master_audit.out" \
        -e "$ERROR_LOG_DIR/master_audit.err" \
        -b y bash src/hpc/run_python_hpc.sh "$PYTHON_EXE" src/python/audit_run_results.py \
            --run-dir "$RUN_BASE_DIR" \
            --raw-dir "$RAW_MASTER_DIR" \
            --output-report "$REPORTS_MASTER_DIR/audit_report_${RUN_NAME}.md" | awk '{print $3}')
    echo "----------------------------------------------------------------------------"
    echo " [9/9] Submitted Automated Post-Run Quality Auditor Job: $audit_job"
    echo "       Will produce audit report on completion: $REPORTS_MASTER_DIR/audit_report_${RUN_NAME}.md"
    echo "----------------------------------------------------------------------------"
fi

echo "============================================================================"
echo " All SGE annotation, filtering, visualization, and audit jobs submitted!"
echo " Check job logs in: $ERROR_LOG_DIR"
echo "============================================================================"
