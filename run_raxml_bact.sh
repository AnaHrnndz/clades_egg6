#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 01-00:00
#SBATCH -e bact_raxml.err
#SBATCH -o bact_raxml.out


export PATH=/home/plaza/soft/standard-RAxML:$PATH

raxmlHPC-PTHREADS -f e -t /home/plaza/projects/eggnog6/ncbi_tree/bifurc_trees/bact_bif_ncbi_tree.nw \
    -s /home/plaza/projects/eggnog6/sp_tree/Bact_sp_tree/cog_all-alg_concat_default-fasttree_default/clean_Bact_COGs_seqs.faa.final_tree.fa \
     -n bact_raxml -m PROTCATLG -T20
