#!/bin/bash

# Logger helper for pipeline execution reporting.
# Sourced by HPC wrapper scripts.

# Ensure LOG_DIR is set
if [ -z "$LOG_DIR" ]; then
    # Fallback to a local logs directory if not set
    LOG_DIR="_log/local_run"
fi
mkdir -p "$LOG_DIR"

REPORT_FILE="$LOG_DIR/pipeline_report.md"
LOCK_FILE="$LOG_DIR/pipeline_report.lock"

# Function to write to report with a file lock to prevent race conditions
write_locked() {
    local content="$1"
    (
        flock -x 200
        echo -e "$content" >> "$REPORT_FILE"
    ) 200>"$LOCK_FILE"
}

# Detect time binary path and return formatting command or empty string if not found
get_time_cmd() {
    local time_log="$1"
    if [ -x "/usr/bin/time" ]; then
        echo "/usr/bin/time -v -o $time_log"
    elif [ -x "/bin/time" ]; then
        echo "/bin/time -v -o $time_log"
    else
        echo ""
    fi
}

# Helper function to run a command timed (using TIME_CMD if configured)
run_command_timed() {
    if [ -n "$TIME_CMD" ]; then
        $TIME_CMD "$@"
    else
        "$@"
    fi
}

# Initialize the report (called at the beginning of the first step or in main.sh)
log_init() {
    local gene="$1"
    local raw_file="$2"
    
    # Check if report already exists, if so reset it
    (
        flock -x 200
        echo "# Pipeline Execution Report for Gene: $gene" > "$REPORT_FILE"
        echo "## General Information" >> "$REPORT_FILE"
        echo "* **Start Date**: $(date '+%Y-%m-%d %H:%M:%S')" >> "$REPORT_FILE"
        echo "* **Host**: $(hostname)" >> "$REPORT_FILE"
        echo "* **Source File**: $raw_file ($(ls -lh "$raw_file" | awk '{print $5}'))" >> "$REPORT_FILE"
        echo -e "\n## Step-by-Step Progress\n" >> "$REPORT_FILE"
    ) 200>"$LOCK_FILE"
}

# Log step start
log_step_start() {
    local step_name="$1"
    local input_file="$2"
    
    # Store start time in epoch seconds for fallback elapsed time calculation
    START_SECS=$(date +%s)
    
    local input_size="N/A"
    if [ -f "$input_file" ]; then
        input_size=$(ls -lh "$input_file" | awk '{print $5}')
    fi
    
    local start_time=$(date '+%Y-%m-%d %H:%M:%S')
    
    local msg="### Step: $step_name
* **Status**: <span style='color: orange;'>Running...</span>
* **Start Time**: $start_time
* **Input File**: \`$(basename "$input_file")\` ($input_size)
* **Allocated Resources**: Memory: ${VMEM:-N/A}, Threads: ${THREADS:-N/A}
"
    write_locked "$msg"
}

# Log step skipped
log_step_skipped() {
    local step_name="$1"
    local output_file="$2"
    
    local output_size="N/A"
    if [ -f "$output_file" ]; then
        output_size=$(ls -lh "$output_file" | awk '{print $5}')
    elif [ -f "${output_file}.gz" ]; then
        output_size=$(ls -lh "${output_file}.gz" | awk '{print $5}')
    fi
    
    local msg="### Step: $step_name
* **Status**: <span style='color: blue;'>Skipped (Reused Existing Output)</span>
* **Output File**: \`$(basename "$output_file")\` ($output_size)
"
    write_locked "$msg"
}

# Log step end (parses time log if available)
log_step_end() {
    local step_name="$1"
    local output_file="$2"
    local exit_code="$3"
    local time_log="$4"
    
    local end_time=$(date '+%Y-%m-%d %H:%M:%S')
    local elapsed="N/A"
    local peak_mem="N/A"
    local output_size="N/A"
    
    # Parse GNU time output if file exists and is not empty
    if [ -f "$time_log" ] && [ -s "$time_log" ]; then
        # Parse maximum resident set size (in KB)
        local max_rss_kb=$(grep "Maximum resident set size" "$time_log" | awk -F': ' '{print $2}' | tr -d ' \t\r\n')
        if [ -n "$max_rss_kb" ]; then
            if [ "$max_rss_kb" -gt 1048576 ]; then
                peak_mem=$(awk "BEGIN {print $max_rss_kb/1048576}")
                peak_mem=$(printf "%.2f GB" "$peak_mem")
            else
                peak_mem=$(awk "BEGIN {print $max_rss_kb/1024}")
                peak_mem=$(printf "%.2f MB" "$peak_mem")
            fi
        fi
        
        # Parse elapsed time
        elapsed=$(grep "Elapsed (wall clock) time" "$time_log" | awk -F': ' '{print $2}' | tr -d ' \t\r\n')
    fi
    
    # Fallback elapsed time if time_log was not generated/valid
    if [ "$elapsed" = "N/A" ] && [ -n "$START_SECS" ]; then
        local end_secs=$(date +%s)
        local diff=$((end_secs - START_SECS))
        if [ "$diff" -ge 3600 ]; then
            elapsed=$(printf "%d:%02d:%02d" $((diff/3600)) $(((diff%3600)/60)) $((diff%60)))
        else
            elapsed=$(printf "%02d:%02d" $((diff/60)) $((diff%60)))
        fi
    fi
    
    # Get output size
    if [ -f "$output_file" ]; then
        output_size=$(ls -lh "$output_file" | awk '{print $5}')
    elif [ -f "${output_file}.gz" ]; then
        output_size=$(ls -lh "${output_file}.gz" | awk '{print $5}')
    fi
    
    local status="<span style='color: green;'>Completed (Success)</span>"
    if [ "$exit_code" -ne 0 ]; then
        status="<span style='color: red;'>Failed (Exit Code $exit_code)</span>"
    fi
    
    local msg="* **End Time**: $end_time
* **Execution Time (wall-clock)**: $elapsed
* **Peak Memory Used**: $peak_mem
* **Output File**: \`$(basename "$output_file")\` ($output_size)
* **Status**: $status
"
    write_locked "$msg"
}

# Compile final log / summary table (called in the final step of the pipeline)
log_final() {
    # We can parse the written report and construct a nice markdown summary table
    (
        flock -x 200
        echo -e "\n## Execution Summary\n" >> "$REPORT_FILE"
        echo "| Step | Status | Start Time | Time | Peak Memory | Output File | Size |" >> "$REPORT_FILE"
        echo "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |" >> "$REPORT_FILE"
        
        # Parse report file to extract step details dynamically
        # Since steps are added in order, we can parse them from the file
        local current_step=""
        local status=""
        local start_time=""
        local elapsed=""
        local peak_mem=""
        local output_file=""
        local size=""
        
        while IFS= read -r line || [ -n "$line" ]; do
            if [[ "$line" =~ ^"### Step: " ]]; then
                # If we had a previous step, write it to the table
                if [ -n "$current_step" ]; then
                    echo "| $current_step | $status | $start_time | $elapsed | $peak_mem | $output_file | $size |" >> "$REPORT_FILE"
                fi
                current_step=$(echo "$line" | sed 's/### Step: //')
                status="N/A"
                start_time="N/A"
                elapsed="N/A"
                peak_mem="N/A"
                output_file="N/A"
                size="N/A"
            elif [[ "$line" =~ ^"* **Status**:" ]]; then
                status=$(echo "$line" | sed 's/\* \*\*Status\*\*: //')
            elif [[ "$line" =~ ^"* **Start Time**:" ]]; then
                start_time=$(echo "$line" | sed 's/\* \*\*Start Time\*\*: //')
            elif [[ "$line" =~ ^"* **Execution Time" ]]; then
                elapsed=$(echo "$line" | sed -E 's/\* \*\*Execution Time.*?\*\*: //')
            elif [[ "$line" =~ ^"* **Peak Memory" ]]; then
                peak_mem=$(echo "$line" | sed 's/\* \*\*Peak Memory Used\*\*: //')
            elif [[ "$line" =~ ^"* **Output File**:" ]]; then
                local raw_out=$(echo "$line" | sed 's/\* \*\*Output File\*\*: //')
                output_file=$(echo "$raw_out" | awk -F'[\` ]' '{print $2}')
                size=$(echo "$raw_out" | awk -F'[()]' '{print $2}')
            fi
        done < "$REPORT_FILE"
        
        # Write the last step
        if [ -n "$current_step" ]; then
            echo "| $current_step | $status | $start_time | $elapsed | $peak_mem | $output_file | $size |" >> "$REPORT_FILE"
        fi
        
        echo -e "\n*Report generated dynamically by pipeline logger on $(date '+%Y-%m-%d %H:%M:%S').*" >> "$REPORT_FILE"
    ) 200>"$LOCK_FILE"
}
