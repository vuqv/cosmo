#!/usr/bin/env bash
# Print the next release version in CalVer YEAR.N form (e.g. 2026.1).
#
#   N = 1 + the highest N among this year's existing "YEAR.N" tags. So the counter
#   auto-increments within a year and auto-RESETS to 1 when the year rolls over
#   (a new year has no "YEAR.*" tags yet). Older "vX.Y" tags are ignored.
#
# Usage:
#   scripts/next_version.sh                 # -> e.g. 2026.1
#   v=$(scripts/next_version.sh) \
#     && git tag -a "$v" -m "release $v" && git push origin "$v"
set -euo pipefail
year="$(date -u +%Y)"
last="$(git tag --list "${year}.*" \
        | sed -n "s/^${year}\.\([0-9]\{1,\}\)$/\1/p" \
        | sort -n | tail -n1)"
printf '%s.%s\n' "$year" "$(( ${last:-0} + 1 ))"
