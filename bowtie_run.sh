#!/bin/bash

input=$1
output=$2
log=$3
name=${input#cutadapt/}
IFS='_' read -r -a array <<< "$name"
genomef=${array[2]}_genome
genome=${array[2]}_bowtie

echo ${genome}
echo ${input}
echo ${output}
(bowtie2 -x genomes/bowtie/${genomef}/${genome} -p 8 --trim-to 43 -q ${input} | samtools view -Sb - > ${output}) 2> ${log}





