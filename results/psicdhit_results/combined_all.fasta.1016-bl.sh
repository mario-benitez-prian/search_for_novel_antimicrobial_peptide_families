#!/bin/sh
#PBS -v PATH
#$ -v PATH


para=$1
cd /home/mario/projects/group_by_family_2/scripts
./../results/psicdhit_results/combined_all.fasta.1016-bl.pl 0 $para &
wait

