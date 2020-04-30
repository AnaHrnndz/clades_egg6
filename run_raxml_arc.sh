#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 01-00:00
#SBATCH -e arc_raxml.err
#SBATCH -o arc_raxml.out


export PATH=/home/plaza/soft/standard-RAxML:$PATH

raxmlHPC-PTHREADS -f e -t /home/plaza/projects/eggnog6/ncbi_tree/bifurc_trees/arc_bif_ncbi_tree.nw -s /home/plaza/projects/eggnog6/sp_tree/Arc_sp_tree/cog_all-alg_concat_default-fasttree_default/Arc2_COG_seqs.faa.final_tree.fa  -n arc_raxml -m PROTCATLG -T20
