#!/usr/bin/env sh
# Yanbin Yin
# 07/21/2015
# hmmscan output parser
# Usage: sh hmmscan-parser.sh hmmscan-output-file
# 1. take hmmer3 --domtblout output and extract necessary columns
# 2. sort on the protein position columns
# 3. remove overlapped/redundant hmm matches and keep the one with the lower e-values
# 4. calculate the covered fraction of hmm &
# apply the E-value cutoff and the covered faction cutoff
cat $1 | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | \
sort -k 3,3 -k 8n -k 9n |
awk '{  if( $5 < 0.001 )  {print}}' | sort | sed 's/\t/,/g'

