#!/usr/bin/bash

MMSEQS=$1
RESULT=$2
THREADS=$3

mkdir -p "${RESULT}"/hhr/
mkdir -p "${RESULT}"/tab/

for NAME in $(ls "${MMSEQS}"/)
do
 hhsearch -i "${MMSEQS}"/"${NAME}" \
 -d /home/niagara/Storage/MetaRus/Common_dir/dbs/hhsuite_dbs/pfama \
 -o "${RESULT}"/hhr/"${NAME}".hhr \
 -blasttab "${RESULT}"/tab/"${NAME}".tab \
 -cpu $THREADS
done

echo 'Finished!'

