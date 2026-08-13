#!/bin/bash
# Universal Genomic Variant Annotation, Filtering, and Visualization Pipeline
# Master Orchestration Script (SGE Cluster Execution)

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

# Customizable Filtering Parameters
MAX_AF=""
MIN_REVEL=""
MIN_ALPHAMISSENSE=""
MIN_SPIP=""
MIN_CADD=""
CONSEQUENCES=""

# Pipeline Skip / Overwrite Flags
skip_variant_converter=false
overwrite_gnomAD=false
skip_merge_gnomad=false
skip_spip_annotation=false
ignore_spip_annotation=false
skip_vep_annotation=false
skip_merge_vep_spip=false
skip_vcf2tsv=false
overwrite_all=false

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
        --skip-spip-annotation)
            skip_spip_annotation=true
            shift 1
            ;;
        --ignore-spip-annotation)
            ignore_spip_annotation=true
            shift 1
            ;;
        --skip-vep-annotation)
            skip_vep_annotation=true
            shift 1
            ;;
        --overwrite-all)
            overwrite_all=true
            shift 1
            ;;
        --sleep-time)
            SLEEP_TIME="$2"
            shift 2
            ;;
        --skip-genes)
            SKIP_GENES="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

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
    echo "  --build <BUILD>            Input genomic assembly build: hg19 or GRCh38 (default: GRCh38)"
    echo "  --output-format <FMT>      Output table format: pq, tsv, xlsx (default: pq)"
    echo "  --max-af <FLOAT>           Filter max gnomAD allele frequency threshold"
    echo "  --min-revel <FLOAT>        Filter min REVEL score threshold"
    echo "  --min-alphamissense <FLOAT> Filter min AlphaMissense score threshold"
    echo "  --min-spip <FLOAT>         Filter min SPiP splice score threshold"
    echo "  --min-cadd <FLOAT>         Filter min CADD phred score threshold"
    echo "  --consequences <LIST>      Comma-separated list of target VEP consequences"
    echo "  --overwrite-all            Re-run all intermediate steps and overwrite results"
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

RAW_FOLDER=$(basename "$RAW_MASTER_DIR")
if [ -z "$RUN_NAME" ]; then
    RUN_NAME="$RAW_FOLDER"
fi

RUN_BASE_DIR="RUNS/${RUN_NAME}"
TMP_MASTER_DIR="${RUN_BASE_DIR}/_tmp"
ANNOTATION_MASTER_DIR="${RUN_BASE_DIR}/annotation"
RESULTS_MASTER_DIR="${RUN_BASE_DIR}/results"
FILTERED_MASTER_DIR="${RUN_BASE_DIR}/filtered"
PLOTS_MASTER_DIR="${RUN_BASE_DIR}/plots"
REPORTS_MASTER_DIR="${RUN_BASE_DIR}/reports"
ERROR_LOG_DIR="${RUN_BASE_DIR}/_log"

GNOMAD_GENES_DIR="annotation/gnomAD_subset"
GENES_BED_FILE="resources/cardio_genes_loc.bed"
GENE_TRANSCRIPT_MAPPING="resources/gene_transcript_mapping.txt"

mkdir -p "$RUN_BASE_DIR" "$TMP_MASTER_DIR" "$ANNOTATION_MASTER_DIR" "$GNOMAD_GENES_DIR" "$RESULTS_MASTER_DIR" "$FILTERED_MASTER_DIR" "$PLOTS_MASTER_DIR" "$REPORTS_MASTER_DIR" "$ERROR_LOG_DIR"

PYTHON_EXE="/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python3"
if [ ! -x "$PYTHON_EXE" ]; then
    PYTHON_EXE="python3"
fi

R_EXE="/data_lab_PGP/shared/utils/conda_envs/datasci/bin/Rscript"
if [ ! -x "$R_EXE" ]; then
    R_EXE="Rscript"
fi

echo "Resolving transcript mappings for target genes..."
$PYTHON_EXE src/python/query_new_transcripts.py "$RAW_MASTER_DIR" --mapping-file "$GENE_TRANSCRIPT_MAPPING" --auto-append || true

hold_all_final_jobs=""

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
    var_conv_job=$(qsub -N "varconv_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 1 \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.varconv.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.varconv.err" \
        -b y bash src/hpc/variant_converter.sh "$input_file" "$vcf_file" --build "$INPUT_BUILD" "$COLUMN_PARAM" | awk '{print $3}')
    echo "  [1/8] Submitted variant_converter job: $var_conv_job"

    # Step 2.1: SPiP Annotation (direct from VCF.gz)
    spip_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annSPiP.vcf"
    spip_job=$(qsub -N "spip_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 4 \
        -hold_jid "$var_conv_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spip.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spip.err" \
        -b y bash src/hpc/annotate_spip_vars.sh "$vcf_gz" "$spip_vcf" | awk '{print $3}')
    echo "  [2/8] Submitted SPiP annotation job: $spip_job"

    # Step 2.2: VEP Annotation (direct from VCF.gz)
    vep_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annVEP.vcf.gz"
    vep_job=$(qsub -N "vep_${gene_name}" -P BIGN -A PGP -l h_vmem=15G -pe smp 4 \
        -hold_jid "$var_conv_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vep.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vep.err" \
        -b y bash src/hpc/annotate_vep_vars.sh "$vcf_gz" "$vep_vcf" | awk '{print $3}')
    echo "  [3/8] Submitted VEP annotation job: $vep_job"

    # Step 2.3: Pangolin Splicing Predictor (direct from VCF.gz)
    pangolin_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annPangolin.vcf.gz"
    pangolin_job=$(qsub -N "pangolin_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 4 \
        -hold_jid "$var_conv_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.pangolin.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.pangolin.err" \
        -b y bash src/hpc/annotate_pangolin_vars.sh "$vcf_gz" "$pangolin_vcf" | awk '{print $3}')
    echo "  [4/8] Submitted Pangolin annotation job: $pangolin_job"

    # Step 2.4: SpliceAI Predictor (Local -D 10000)
    spliceai_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annSpliceAI.vcf.gz"
    spliceai_job=$(qsub -N "spliceai_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 4 \
        -hold_jid "$var_conv_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spliceai.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.spliceai.err" \
        -b y bash src/hpc/annotate_spliceai_vars.sh "$vcf_gz" "$spliceai_vcf" | awk '{print $3}')
    echo "  [5/8] Submitted SpliceAI (-D 10000) annotation job: $spliceai_job"

    # Step 2.5: Branchpoint Predictor Pair (Branchpointer + LaBranchoR)
    branchpoint_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annBranchpoint.vcf.gz"
    branchpoint_job=$(qsub -N "branchpoint_${gene_name}" -P BIGN -A PGP -l h_vmem=10G -pe smp 2 \
        -hold_jid "$var_conv_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.branchpoint.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.branchpoint.err" \
        -b y bash src/hpc/annotate_branchpointer_vars.sh "$vcf_gz" "$branchpoint_vcf" | awk '{print $3}')
    echo "  [6/8] Submitted Branchpoint annotation job: $branchpoint_job"

    # Step 3: Merge all Predictor Annotations into Final Annotated VCF
    final_annotated_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.annotated.vcf.gz"
    large_sv_vcf="$ANNOTATION_MASTER_DIR/${gene_name}.large_svs.vcf.gz"
    merge_ann_job=$(qsub -N "merge_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
        -hold_jid "$vep_job,$spip_job,$pangolin_job,$spliceai_job,$branchpoint_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.merge.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.merge.err" \
        -b y bash src/hpc/merge_vep_spip.sh "$spip_vcf.gz" "$vep_vcf" "$final_annotated_vcf" "$large_sv_vcf" "$pangolin_vcf" "$spliceai_vcf" "$branchpoint_vcf" | awk '{print $3}')
    echo "  [7/8] Submitted multi-predictor merge job: $merge_ann_job"

    # Step 4: Direct VCF to Final Clean Table (Parquet/Excel/TSV)
    vcf2parsed_job=$(qsub -N "vcf2parsed_${gene_name}" -P BIGN -A PGP -l h_vmem=80G -pe smp 1 \
        -hold_jid "$merge_ann_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vcf2parsed.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.vcf2parsed.err" \
        -b y bash src/hpc/vcf2parsed.sh "$final_annotated_vcf" "$final_output_file" "$Transcript" | awk '{print $3}')
    echo "  [8/8] Submitted direct vcf2parsed job: $vcf2parsed_job"

    # Step 5: Downstream Filtering & Graph Generation
    filter_flags="--input $final_output_file --output-dir $FILTERED_MASTER_DIR"
    [ -n "$MAX_AF" ] && filter_flags="$filter_flags --max-af $MAX_AF"
    [ -n "$MIN_REVEL" ] && filter_flags="$filter_flags --min-revel $MIN_REVEL"
    [ -n "$MIN_ALPHAMISSENSE" ] && filter_flags="$filter_flags --min-alphamissense $MIN_ALPHAMISSENSE"
    [ -n "$MIN_SPIP" ] && filter_flags="$filter_flags --min-spip $MIN_SPIP"
    [ -n "$MIN_CADD" ] && filter_flags="$filter_flags --min-cadd $MIN_CADD"
    [ -n "$CONSEQUENCES" ] && filter_flags="$filter_flags --consequences \"$CONSEQUENCES\""

    filter_job=$(qsub -N "filter_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
        -hold_jid "$vcf2parsed_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.filter.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.filter.err" \
        -b y $PYTHON_EXE src/python/filter_and_summarize.py $filter_flags | awk '{print $3}')

    plot_job=$(qsub -N "plot_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
        -hold_jid "$filter_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.plot.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.plot.err" \
        -b y $R_EXE src/R/plot_annotation_results.R --input "$final_output_file" --output-dir "$PLOTS_MASTER_DIR" | awk '{print $3}')

    report_job=$(qsub -N "report_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
        -hold_jid "$filter_job" \
        -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.report.out" \
        -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.report.err" \
        -b y $PYTHON_EXE src/python/generate_interactive_report.py --input "$final_output_file" --output "$REPORTS_MASTER_DIR/${gene_name}_interactive_dashboard.html" | awk '{print $3}')

    sleep $SLEEP_TIME
done

echo "============================================================================"
echo " All SGE annotation, filtering, and visualization jobs submitted!"
echo " Check job logs in: $ERROR_LOG_DIR"
echo "============================================================================"
