#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: scripts/run_batch.sh FASTQ_DIR RUN_ID

Environment overrides:
  THREADS=8
  UNIT=batch
  JOBS=1                 # number of samples to run at once
  OUT_ROOT=analysis/nips_viral
  RESOURCES=$PWD/resources
  IMAGE=localhost/nips_viral_screen:v1.0.0
  PLATFORM=linux/amd64
  SAMPLE_REGEX=          # optional regex matched against parsed sample names
  RESUME_COMPLETED=1     # skip samples with an existing viral JSON
  FORCE=0                # rerun completed samples when set to 1

Input filenames supported:
  SAMPLE.R1.fastq.gz / SAMPLE.R2.fastq.gz
  SAMPLE_S42_R1_001.fastq.gz / SAMPLE_S42_R2_001.fastq.gz
  SAMPLE_S42_L001_R1_001.fastq.gz / SAMPLE_S42_L001_R2_001.fastq.gz
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -ne 2 ]]; then
  usage
  exit 1
fi

fastq_dir=$(cd "$1" && pwd)
run_id=$2

threads=${THREADS:-8}
unit=${UNIT:-batch}
jobs=${JOBS:-1}
out_root=${OUT_ROOT:-analysis/nips_viral}
resources=${RESOURCES:-$PWD/resources}
image=${IMAGE:-localhost/nips_viral_screen:v1.0.0}
platform=${PLATFORM:-linux/amd64}
sample_regex=${SAMPLE_REGEX:-}
resume_completed=${RESUME_COMPLETED:-1}
force=${FORCE:-0}

if ! [[ "$jobs" =~ ^[0-9]+$ ]] || [[ "$jobs" -lt 1 ]]; then
  echo "JOBS must be a positive integer" >&2
  exit 2
fi
if [[ "$resume_completed" != 0 && "$resume_completed" != 1 ]]; then
  echo "RESUME_COMPLETED must be 0 or 1" >&2
  exit 2
fi
if [[ "$force" != 0 && "$force" != 1 ]]; then
  echo "FORCE must be 0 or 1" >&2
  exit 2
fi

if [[ ! -d "$resources" ]]; then
  echo "Missing resources directory: $resources" >&2
  exit 2
fi
resources=$(cd "$resources" && pwd)

if [[ ! -d "$resources/Homo_sapiens/hg19" ]]; then
  echo "Missing resources/Homo_sapiens/hg19 under: $resources" >&2
  exit 2
fi

run_root="$out_root/$run_id"
mkdir -p "$run_root"
run_root=$(cd "$run_root" && pwd)
samples_root="$run_root/samples"
work_root="$run_root/work"
logs_root="$run_root/logs"
json_root="$run_root/summary/viral_json"
reports_root="$run_root/reports"
mkdir -p "$samples_root" "$work_root" "$logs_root" "$json_root" "$reports_root"
r1_list="$run_root/r1_files.txt"
manifest="$run_root/manifest.tsv"

sample_name_for() {
  local base=$1
  if [[ $base =~ ^(.+)_S[0-9]+(_L[0-9]{3})?_R1_001\.f(ast)?q\.gz$ ]]; then
    printf '%s\n' "${BASH_REMATCH[1]}"
  elif [[ $base =~ ^(.+)\.R1\.f(ast)?q\.gz$ ]]; then
    printf '%s\n' "${BASH_REMATCH[1]}"
  else
    printf '%s\n' "${base%%_R1*}"
  fi
}

r2_for() {
  local r1=$1
  local r2=${r1/_R1_001./_R2_001.}
  r2=${r2/.R1./.R2.}
  printf '%s\n' "$r2"
}

stage_file() {
  local src=$1
  local dest=$2
  local tmp="$dest.tmp.$$"
  local src_size
  local dest_size

  if [[ -e "$dest" ]]; then
    src_size=$(wc -c < "$src" | tr -d ' ')
    dest_size=$(wc -c < "$dest" | tr -d ' ')
    if [[ "$src_size" == "$dest_size" ]]; then
      return 0
    fi
    echo "Replacing incomplete staged file: $dest" >&2
    rm -f "$dest"
  fi

  rm -f "$tmp"
  cp "$src" "$tmp"
  mv "$tmp" "$dest"
}

copy_run_artifacts() {
  local sample=$1
  local sample_dir=$2
  local json_file="$sample_dir/nips_viral/v1.0.0/$sample.viral.json"

  [[ -f "$json_file" ]] && cp "$json_file" "$json_root/"
  [[ -f "$sample_dir/nextflow_trace.txt" ]] && cp "$sample_dir/nextflow_trace.txt" "$reports_root/$sample.trace.txt"
  [[ -f "$sample_dir/nextflow_report.html" ]] && cp "$sample_dir/nextflow_report.html" "$reports_root/$sample.report.html"
  [[ -f "$sample_dir/nextflow_timeline.html" ]] && cp "$sample_dir/nextflow_timeline.html" "$reports_root/$sample.timeline.html"
}

sample_completed() {
  local sample=$1
  local sample_dir=$2
  local json_file="$sample_dir/nips_viral/v1.0.0/$sample.viral.json"

  [[ -s "$json_file" ]]
}

find "$fastq_dir" -maxdepth 1 -type f \
  \( -name '*_R1_001.fastq.gz' -o -name '*.R1.fastq.gz' -o -name '*_R1_001.fq.gz' -o -name '*.R1.fq.gz' \) \
  | sort > "$r1_list"

if [[ ! -s "$r1_list" ]]; then
  echo "No R1 FASTQ files found in $fastq_dir" >&2
  exit 1
fi

printf 'sample\tsource_r1\tsource_r2\tstaged_r1\tstaged_r2\n' > "$manifest"

run_sample() {
  local sample=$1
  local sample_dir=$2
  local sample_work_dir=$3
  local log_file=$4
  local json_file

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running $sample"
  docker run --platform "$platform" --rm \
    -v "$sample_dir:/tmp/sample:rw" \
    -v "$resources:/tmp/ref:ro" \
    -v "$sample_work_dir:/tmp/work:rw" \
    "$image" \
    --resources /tmp/ref \
    --sample_path /tmp/sample \
    --run "$run_id" \
    --samplename "$sample" \
    --unit "$unit" \
    --threads "$threads" \
    -work-dir /tmp/work \
    -with-trace "/tmp/sample/nextflow_trace.txt" \
    -with-report "/tmp/sample/nextflow_report.html" \
    -with-timeline "/tmp/sample/nextflow_timeline.html" 2>&1 | tee "$log_file"

  json_file="$sample_dir/nips_viral/v1.0.0/$sample.viral.json"
  if [[ -f "$json_file" ]]; then
    copy_run_artifacts "$sample" "$sample_dir"
  else
    echo "Expected JSON not found after successful run: $json_file" >&2
    return 1
  fi
}

pids=()
labels=()
failures=0
selected=0

wait_first() {
  local pid=${pids[0]}
  local label=${labels[0]}
  local status=0

  wait "$pid" || status=$?
  if [[ "$status" -ne 0 ]]; then
    echo "Sample failed: $label" >&2
    failures=1
  fi

  pids=("${pids[@]:1}")
  labels=("${labels[@]:1}")
  return "$status"
}

wait_all() {
  while [[ ${#pids[@]} -gt 0 ]]; do
    wait_first || true
  done
}

while IFS= read -r r1; do
  base=$(basename "$r1")
  sample=$(sample_name_for "$base")
  if [[ -n "$sample_regex" && ! "$sample" =~ $sample_regex ]]; then
    continue
  fi
  selected=$((selected + 1))

  r2=$(r2_for "$r1")
  if [[ ! -f "$r2" ]]; then
    echo "Missing R2 pair for $r1; expected $r2" >&2
    exit 1
  fi

  sample_dir="$samples_root/$sample"
  sample_fastq_dir="$sample_dir/fastq"
  sample_work_dir="$work_root/$sample"
  log_file="$logs_root/$sample.log"

  mkdir -p "$sample_fastq_dir" "$sample_work_dir"
  stage_file "$r1" "$sample_fastq_dir/$sample.R1.fastq.gz"
  stage_file "$r2" "$sample_fastq_dir/$sample.R2.fastq.gz"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$sample" \
    "$r1" \
    "$r2" \
    "$sample_fastq_dir/$sample.R1.fastq.gz" \
    "$sample_fastq_dir/$sample.R2.fastq.gz" >> "$manifest"

  if [[ "$force" -eq 0 && "$resume_completed" -eq 1 ]] && sample_completed "$sample" "$sample_dir"; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping completed sample $sample"
    copy_run_artifacts "$sample" "$sample_dir"
    continue
  fi

  run_sample "$sample" "$sample_dir" "$sample_work_dir" "$log_file" &
  pids+=("$!")
  labels+=("$sample")

  if [[ ${#pids[@]} -ge "$jobs" ]]; then
    wait_first || true
  fi

  if [[ "$failures" -ne 0 ]]; then
    break
  fi
done < "$r1_list"

wait_all

if [[ "$selected" -eq 0 ]]; then
  echo "No samples matched SAMPLE_REGEX=$sample_regex" >&2
  exit 1
fi

if [[ "$failures" -ne 0 ]]; then
  echo "One or more samples failed. Check logs in: $logs_root" >&2
  exit 1
fi

echo "Done. Per-sample outputs: $samples_root"
echo "JSON summary copies: $json_root"
echo "Nextflow reports: $reports_root"
echo "Manifest: $manifest"
