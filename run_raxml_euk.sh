#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 01-00:00
#SBATCH -e euk_raxml.err
#SBATCH -o euk_raxml.out


export PATH=/home/plaza/soft/standard-RAxML:$PATH

raxmlHPC-PTHREADS -f e -t /home/plaza/projects/eggnog6/ncbi_tree/bifurc_trees/euk_bif_ncbi_tree.nw \
-s /home/plaza/projects/eggnog6/sp_tree/Euk_sp_tree/cog_all-alg_concat_default-fasttree_full/proteome_seqs.faa.final_tree.fa -n euk_raxml -m PROTCATLG -T20
