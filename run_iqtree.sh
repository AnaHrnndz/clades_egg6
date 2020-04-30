#!/bin/bash
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 5-00:00
#SBATCH -e iqtree.err
#SBATCH -o iqtree.out


export PATH=/home/plaza/soft/iqtree-2.0-rc1-Linux/bin:$PATH

iqtree -s /home/plaza/projects/eggnog6/sp_tree_marker_genes/basic_sptree/cog_all-alg_concat_default-fasttree_default/prune_COGs_seqs.fasta.final_tree.fa -m LG \
-g /home/plaza/projects/eggnog6/guide_tree_ncbi/post_process_freeze12_representative_sp_tree_oldDB.nw -nt 40 -fast -sprrad 1 -ninit 10
