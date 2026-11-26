#!/usr/bin/env bash
# Run the benchmark at each scale in its own process and record process peak RSS.
# /usr/bin/time wraps Rscript inside the pixi env so it measures R itself, not pixi.
# Usage (from package root): bash bench/run_all.sh [label] [scales...]
#   bash bench/run_all.sh baseline S M L
set -euo pipefail

label="${1:-$(date +%Y-%m-%d)}"
shift || true
if [ "$#" -eq 0 ]; then scales=(S M L); else scales=("$@"); fi

PIXI="${PIXI:-$HOME/.pixi/bin/pixi}"
mkdir -p bench/results

for s in "${scales[@]}"; do
  echo "=== scale $s ==="
  log="bench/results/${label}_${s}.time.log"
  # macOS /usr/bin/time -l reports "maximum resident set size" in bytes;
  # GNU time -v (Linux) reports "Maximum resident set size (kbytes)".
  if [[ "$(uname)" == "Darwin" ]]; then
    "$PIXI" run /usr/bin/time -l Rscript bench/run_benchmarks.R "$s" "$label" 2> >(tee "$log" >&2)
    rss_mb=$(awk '/maximum resident set size/ {printf "%.0f", $1/1048576}' "$log")
  else
    "$PIXI" run /usr/bin/time -v Rscript bench/run_benchmarks.R "$s" "$label" 2> >(tee "$log" >&2)
    rss_mb=$(awk -F: '/Maximum resident set size/ {printf "%.0f", $2/1024}' "$log")
  fi
  echo "peak RSS for scale $s: ${rss_mb} MB"
  echo "${label},${s},${rss_mb}" >> bench/results/peak_rss.csv
done
