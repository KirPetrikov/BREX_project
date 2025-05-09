#!/usr/bin/bash

MMSEQS=$1
RESULT=$2
THREADS=$3

mkdir -p "${RESULT}"

for NAME in $(ls "${MMSEQS}"/)
do
 hhblits -i "${MMSEQS}"/${NAME} \
 -o "${RESULT}"/${NAME}.hhr \
 -d /home/niagara/Storage/MetaRus/Common_dir/dbs/hhsuite_dbs/pdb70 \
 -cpu $THREADS -n 1
done

