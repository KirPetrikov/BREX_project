#!/bin/bash
#v0.1
# $1 - path to input fasta-file
# $2 - path to save results
# $3 - N threads

FASTA=$1
RESULT=$2
THREADS=$3

mkdir -p "${RESULT}"/db

mmseqs createdb "${FASTA}" "${RESULT}"/db/db

mkdir "${RESULT}"/clust

mmseqs cluster "${RESULT}"/db/db "${RESULT}"/clust/clust "${RESULT}"/mmseqs_tmp --min-seq-id 0.30 -c 0.8 --cov-mode 0 --threads $THREADS

mmseqs createtsv "${RESULT}"/db/db "${RESULT}"/db/db "${RESULT}"/clust/clust "${RESULT}"/mmseqs_clusters.tsv --threads $THREADS

echo '*** Create cluster table and summary histogram ***'

python3 mmseqs_clusters_table.py -i "${RESULT}"/mmseqs_clusters.tsv -o "${RESULT}"

echo '*** Finished ***'

