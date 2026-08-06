#!/bin/bash
# Show GPML controller and child SLURM jobs, including Snakemake rule details.
set -euo pipefail

user="${1:-$USER}"
command -v squeue >/dev/null 2>&1 || { echo "ERROR: squeue was not found in PATH." >&2; exit 1; }

echo "Active SLURM jobs for ${user} (DETAIL contains Snakemake rule/wildcards when available)"
squeue -u "$user" -h -o "%i|%t|%C|%m|%M|%l|%j|%R|%k" | awk -F '|' '
function gib(value, number, unit) { number=value; sub(/[A-Za-z]+$/, "", number); unit=value; sub(/^[0-9.]+/, "", unit); if (unit == "K") return sprintf("%.1f", number/1024/1024); if (unit == "M") return sprintf("%.1f", number/1024); if (unit == "G" || unit == "") return sprintf("%.1f", number); if (unit == "T") return sprintf("%.1f", number*1024); return value }
BEGIN { printf "%-12s %-2s %5s %9s %10s %11s %-30s %-18s %s\\n", "JOBID", "ST", "CPUS", "MEM_GIB", "TIME", "TIME_LIMIT", "NAME", "NODE/REASON", "DETAIL" }
{ detail=($9 == "" || $9 == "(null)") ? "-" : $9; printf "%-12s %-2s %5s %9s %10s %11s %-30s %-18s %s\\n", $1,$2,$3,gib($4),$5,$6,substr($7,1,30),substr($8,1,18),detail }'
