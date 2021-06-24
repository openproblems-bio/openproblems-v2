#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'modality_alignment|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

TARGET=target/docker/modality_alignment
OUTPUT=output_bash/modality_alignment

mkdir -p $OUTPUT/datasets
mkdir -p $OUTPUT/methods
mkdir -p $OUTPUT/metrics

# generate datasets
if [ ! -f "$OUTPUT/datasets/citeseq_cbmc.h5ad" ]; then
  tmp1=`tempfile`
  tmp2=`tempfile`
  wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz' -O "$tmp1"
  wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz' -O "$tmp2"
  "$TARGET/datasets/scprep_csv/scprep_csv" \
    --input1 "$tmp1" \
    --input2 "$tmp2" \
    --output "$OUTPUT/datasets/citeseq_cbmc.h5ad"
  rm "$tmp1" "$tmp2"
fi

# run all methods on all datasets
for meth in `ls "$TARGET/methods"`; do
  for dat in `ls "$OUTPUT/datasets"`; do
    dat_id="${dat%.*}"
    input_h5ad="$OUTPUT/datasets/$dat_id.h5ad"
    output_h5ad="$OUTPUT/methods/${dat_id}_$meth.h5ad"
    if [ ! -f "$output_h5ad" ]; then
      echo "> $TARGET/methods/$meth/$meth -i $input_h5ad -o $output_h5ad"
      "$TARGET/methods/$meth/$meth" -i "$input_h5ad" -o "$output_h5ad"
    fi
  done
done

# run all metrics on all outputs
for met in `ls "$TARGET/metrics"`; do
  for outp in `ls "$OUTPUT/methods"`; do
    out_id="${outp%.*}"
    input_h5ad="$OUTPUT/methods/$out_id.h5ad"
    output_h5ad="$OUTPUT/metrics/${out_id}_$met.h5ad"
    if [ ! -f "$output_h5ad" ]; then
      echo "> $TARGET/metrics/$met/$met" -i "$input_h5ad" -o "$output_h5ad"
      "$TARGET/metrics/$met/$met" -i "$input_h5ad" -o "$output_h5ad"
    fi
  done
done

# concatenate all scores into one tsv
INPUTS=$(ls -1 "$OUTPUT/metrics" | sed "s#.*#-i '$OUTPUT/metrics/&'#" | tr '\n' ' ')
eval "$TARGET/../common/extract_scores/extract_scores" $INPUTS -o "$OUTPUT/scores.tsv"
