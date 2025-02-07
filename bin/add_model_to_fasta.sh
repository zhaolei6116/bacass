#!/bin/bash
MODEL=$1
FASTA=$2
sed -i "/>/ s/$/ basecall_model=$MODEL/" $FASTA
