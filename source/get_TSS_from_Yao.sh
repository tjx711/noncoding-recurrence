#!/bin/bash

# ../src/bash/get_TSS_metadata_from_UCSC.sh csnps.rslist.txt csnps.tss_metadata.txt
# ../src/bash/get_TSS_metadata_from_UCSC.sh rsnps.rslist.txt rsnps.tss_metadata.txt

# Process Cmd Line Args
if [ $# -lt 2 ]
then
  echo "Usage: $0 rslist.txt tss_metadata.txt"
  exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2
INPUT_CHUNK_SIZE=1000

# Check if ${INPUT_FILE} exists

if [ ! -f "${INPUT_FILE}" ]
then
	echo "Input file ${INPUT_FILE} does not exists."
	exit -1
fi

# If ${OUTPUT_FILE} exists, delete it,
# because we are going to APPEND the results from queries to it.

if [ -f "${OUTPUT_FILE}" ]
then
	rm ${OUTPUT_FILE}
	echo "A previous output file ${OUTPUT_FILE} is deleted."
fi

# Query

function queryRemoteUcscForTss {
    flag_skip_column=""

	# if file exists, do not append header
	if [ -f $2 ]
	then
		flag_skip_column="-N"
	fi

	# ${INPUT_FILE}
		# rs1
		# rs2
		# ...
	# cat ${INPUT_FILE} | xargs -n1000
		# rs1 rs2 ...
	# SNPS=\'${1// /\',\'}\'
		# where $1 is from `cat ${INPUT_FILE} | xargs -n1000`
		# SNPS="'rs1','rs2',..."
		# Add single quotes and commas for `WHERE r.name IN(...)`

	# Replace every inbetween space by `','`; add `'` to both head and tail
    SNPS=\'${1// /\',\'}\'

    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A ${flag_skip_column} -e "\
        SET sql_mode = 'NO_UNSIGNED_SUBTRACTION';
        SELECT
            s.name, s.chrom, s.chromStart, s.chromEnd,
            g.name as tssGene, g.txStart, g.txEnd, g.strand,
        CASE g.strand
            WHEN '+' THEN g.txStart - s.chromStart
            WHEN '-' THEN s.chromStart - g.txEnd   # Steve edited so that 'upstream' is always a positive distance, 2015.10.08
            # WHEN '-' THEN g.txEnd - s.chromStart  # Satpreet's original code 2015.08.31
        END as tssDistance
        FROM
            snp142 s
        LEFT OUTER JOIN
            ensGene g
        ON
            g.bin = s.bin # Speeds up JOINs
            AND g.chrom = s.chrom
        WHERE
            s.name IN  ( ${SNPS} )
        ORDER BY
            name, abs(tssDistance)
    " hg19 >> $2
}

# export makes the variable available to sub-processes
	# http://stackoverflow.com/a/1158231
export -f queryRemoteUcscForTss

cat ${INPUT_FILE} | xargs -n${INPUT_CHUNK_SIZE} | xargs -I {} bash -c 'queryRemoteUcscForTss "$@"' _ {} ${OUTPUT_FILE}

python ../python/getClosestTSS.py ${OUTPUT_FILE} # transforms file in-place
