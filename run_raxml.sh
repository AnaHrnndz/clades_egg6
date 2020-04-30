#!/bin/bash
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 01-00:00
#SBATCH -e raxml.err
#SBATCH -o raxml.out


export PATH=/home/plaza/soft/standard-RAxML:$PATH

raxmlHPC-PTHREADS -f e -t /home/plaza/projects/eggnog6/guide_tree_ncbi/post_process_freeze12_representative_sp_tree_oldDB.nw -s /home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/prune_COGs_seqs.fasta.final_tree.fa -n branch -m PROTCATLG -T20
