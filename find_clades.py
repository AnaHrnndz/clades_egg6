import sys
import os
import numpy
import itertools

from ete3 import Tree, NCBITaxa
from ete3.utils import print_table
from argparse import ArgumentParser

LEVELS_FILE = "/home/plaza/projects/eggnog6/2_bact_newXNogs_egg6.tsv"
TAXONOMY_DB = "/home/plaza/projects/eggnog6/ncbi_taxa_db/taxa.sqlite"

def distance_matrix(t):
    dist2root = {}
    current_dist = 0
    node2node = {}
    content = t.get_cached_content()
    for post, node in t.iter_prepostorder():
        if node.dist < 0:
            node.dist = 0.0
        if post:
            for ch1, ch2 in itertools.combinations(node.children, 2):
                content1 = content[ch1]
                content2 = content[ch2]                    
                for leaf1, leaf2 in itertools.product(content1, content2):
                    d1 = dist2root[leaf1] - current_dist
                    d2 = dist2root[leaf2] - current_dist
                    node2node[(leaf1, leaf2)] = d1 + d2
                    node2node[(node, leaf1)] = d1
                    node2node[(node, leaf2)] = d2                
                    
            current_dist -= node.dist                
        else:
            dist2root[node] = current_dist + node.dist
            if not node.is_leaf():
                current_dist += node.dist
    return content, node2node

def scan_level(level_node, content, node2node_dist, thr, min_clade_factor, dev_thr, forced_clades, size_penalty):
    indent = 2

    # Calculate some values about the level node
    leaf_to_level_distances = [node2node_dist[(level_node, leaf)] for leaf in content[level_node]]
    level_leaf_distance_mean = numpy.mean(leaf_to_level_distances)
    level_leaf_distance_median = numpy.median(leaf_to_level_distances)
    level_leaf_distance_max = numpy.max(leaf_to_level_distances)
    level_leaf_distance_deviation = numpy.std(leaf_to_level_distances)
    nleaves = len(leaf_to_level_distances)
    allow_single_species = False
    
    # This is the threshold function    
    def is_leaf_based_on_distance_to_level(_node):
        ''' Should return True if a node is considered a clade under current level'''
        if not allow_single_species and len(content[_node]) == len(_node.children):
            return True
        
        if _node is level_node:
            return False
        
        if _node.is_leaf():
            return True

        if _node in forced_clades:
            return True
            
        # This is the distance of the potential clade to the parent level node
        _dist_to_level = _node.get_distance(level_node)
        
        #if node is sufficiently distant from the level node, reported as a clade node
        if  _dist_to_level >= level_leaf_distance_max * (1-relaxation):
            return True
        else:
            return False
        
    thr_function = is_leaf_based_on_distance_to_level
    level_sci_name = sci_names[int(level_node.name)]
    print ()
    print ("** Level name:                        ", level_sci_name, level_node.name)
    print ("   Mean leaf distance to level node:  ", level_leaf_distance_mean)
    print ("   Median leaf distance to level node:", level_leaf_distance_median)
    print ("   Std. leaf distance to level node:  ", level_leaf_distance_deviation)
    print ("   Max leaf distance to level node:   ", level_leaf_distance_max)
    print ("   Level Size:                        ", len(content[level_node]))
    relaxation = 1.0
    while True:        
        selected_clades = {}
        for nclades, node in enumerate(level_node.iter_leaves(is_leaf_fn=thr_function)):
            D = node.get_distance(level_node)
            if node.children:
                leaf_to_clade_distances = [node2node_dist[(node, _leaf)] for _leaf in content[node]]
                clade_leaf_distance_mean = numpy.mean(leaf_to_clade_distances)
                clade_leaf_distance_median = numpy.median(leaf_to_clade_distances)
                clade_leaf_distance_deviation  = numpy.std(leaf_to_clade_distances)
            else:
                clade_leaf_distance_mean = node.dist
                clade_leaf_distance_median = node.dist
                clade_leaf_distance_deviation  = node.dist
                
            leaf_names = [leaf.name for leaf in content[node]] if node.children else [node.name]            
            selected_clades[node.name] = [D, len(content[node]), clade_leaf_distance_mean, clade_leaf_distance_median, clade_leaf_distance_deviation, leaf_names, node]
            
        #min_clade_factor_corrected = min_clade_factor * (1 - ((nleaves / LARGEST_CLADE_SIZE) * (1-size_penalty))
        expected_size = nleaves * min_clade_factor * (1 - ((nleaves / LARGEST_CLADE_SIZE) * size_penalty))
        min_size = max(3, round(expected_size))
               
        if len(selected_clades) < min_size:
            # relax the splitting thr to get more clades in the next round
            relaxation -= 0.01
            relaxation = max(relaxation, 0.0)
            if relaxation == 0:
                allow_single_species = True
            #print selected_clades
            #print_clade_tree(selected_clades)
            print ("   Relaxing thr to get at least %s clades: %s %s dist-cutoff:%s %s" %(min_size, level_sci_name, relaxation, level_leaf_distance_max * (1-relaxation), len(selected_clades)))
        # elif len(selected_clades) == len(content[level_node]):
        #     relaxation += 0.01
        #     print "Number of clades equal to size, increasing thr", level_sci_name, relaxation
        else:
            break
    #print "** SELECTION MIN SIZE", min_size, nleaves * min_clade_factor , expected_size, LARGEST_CLADE_SIZE       
    return selected_clades
            
def print_clade_tree(selected_clades, level=None, level_name=None):
    print ('selected_clades', selected_clades.keys())
    clade_tree = ncbi.get_topology(selected_clades.keys())
    tax2name, tax2lin, tax2rank = ncbi.annotate_tree(clade_tree)
    _, syns = ncbi._translate_merged(selected_clades)
    if syns: 
        new2old = dict([(v,k) for k,v in syns.items()])
        for n in clade_tree:
            n.name = str(new2old.get(int(n.name), n.name))
    
    table = []
    header = "taxid name size dist_to_level_node mean_dists meadian_dists dev_dists".split()
    total_size = 0
    for clade, (D, size, dmean, dmedian, ddev, leaves, clade_node) in selected_clades.items():
        table.append([clade, tax2name.get(int(clade), clade), size, D, dmean, dmedian, ddev])
        total_size += size
        
    table_title = "%s (%d Leaves %d selected clades)" %(level_name, total_size, len(selected_clades))
    for node in clade_tree:
        node.D = selected_clades[node.name][0]
        node.S = selected_clades[node.name][1]
        node.C = selected_clades[node.name][2]
        
    print (clade_tree.get_ascii(attributes=["sci_name", "D", "S", "C"]))
    print_table(table, max_col_width=60, title=table_title, header=header)    
    print ("** SELECTION %s %s: %d/%d" %(level, level_name, len(selected_clades), total_size))

    if args.output:
        filename = os.path.join(args.output, "%s.%s.clades" %(level, level_name.replace('/', '_').replace(' ', "_")))        
        with open(filename, "w") as OUT:
            for clade, (D, size, dmean, dmedian, ddev, leaves, clade_node) in selected_clades.items():
                print ('\t'.join(map(str, [clade, size, tax2name.get(int(clade), clade), ','.join(leaves)])), file = OUT,)
                              
def directory(name):
    if os.path.isdir(name):
        return name
    else:
        raise ValueError("Path to directory expected")
                
if __name__ == "__main__":
    parser = ArgumentParser("select clades")

    parser.add_argument("-t", dest="guide_tree", help="guiding tree", required=True)

    parser.add_argument("-m", dest="min_clades", default=0.25, type=float, help='node size factor from 0 to 1 used to establish the min number of clades that should be produced for each level. A factor=0.33 would generate ~33 clades for a level with 100 species. Default=0.25')
    parser.add_argument("-p", dest="size_penalty", default=0.0, type=float, help='size penality applied to the size factor. Default=0.0' )
    parser.add_argument("-l", dest="start_level", help="starting level point (taxid). All sub levels will also be processed")   
    parser.add_argument("--output", dest="output", help="write per-level clades information into files under the provided path.", type=directory)
    parser.add_argument("--ignore_distances", dest="ignore_distances", help="do not use real branch lengths, only node depth")
    
    # parser.add_argument("-s", dest="dev_thr", help="Call clades based on internal branch deviation within clade topology.  1=splitting is easy (more clades per levek), 0=splitting is difficult(less clades per level)", default=0.1, type=float)
    # parser.add_argument("-d", dest="distance_factor", help="Call clades based in their distance to the level node. 1=splitting is easy (more clades per levek), 0=splitting is difficult(less clades per level)", default=0.4, type=float)
    # parser.add_argument("-f", dest="level_file", help="file", default="")
   
    args = parser.parse_args()

    ncbi= NCBITaxa(TAXONOMY_DB)
    
    # Load tree with internal node names
    t = Tree(args.guide_tree, format=1)
    ncbi.annotate_tree(t, taxid_attr='name')
    
    for n in t.traverse():
        n.name = str(n.taxid)
        if n.dist < 0:
            n.dist = 0
        if args.ignore_distances: 
            n.dist = 1            
    
    # load predefined levels and indentify nodes in target guiding tree
    print ('Processing levels file', LEVELS_FILE)
    levels = set()
    level2species = {}
    expected_taxa = set()
    level_names = {}
    for line in open(LEVELS_FILE):
        line = line.strip()
        if line.startswith("#"):
            continue
        #level_taxid, iname, lname, ntaxa, level_taxa_raw = line.split('\t')
        if len(line.split('\t')) >3:
            level_taxid,  lname, ntaxa, level_taxa_raw = line.split('\t')
        else:
            continue
        
        level_taxa = map(str.strip, level_taxa_raw.split())
        level_taxa_cp = map(str.strip, level_taxa_raw.split())
        
        level_names[level_taxid] = lname
        
        #expected_taxa.update(level_taxa)
        expected_taxa.update(level_taxa_cp)
        
        levels.add(level_taxid.strip())
        
        level2species[level_taxid.strip()] = set(level_taxa)



    print (len(expected_taxa - set(t.get_leaf_names())), "missing species")
    
    #Fix names in tree to match level names
    root_node = t
    levels_found = 0
    for lv, taxa in level2species.items():
        #print(lv, len(list(taxa)))
        if len(taxa) > 1:
            n = t.get_common_ancestor(taxa)
        else:
            continue
            #n = t & list(taxa)[0]
            
        if len(n) == len(taxa):
            levels_found += 1
            n.name = lv
            if args.start_level and lv == args.start_level:
                root_node = n
        else:
            print ("ERROR", lv, level_names[lv], len(taxa), len(n), set(taxa) ^ set(n.get_leaf_names()))
            
        sys.stdout.flush()
    print ("** %d target levels found in tree" %levels_found)

    # Use NCBI names for better debugging
    print ("Loading NCBI data..." )
    complete_ref = ncbi.get_topology(t.get_leaf_names(), intermediate_nodes=True)    
    try:
        sci_names = ncbi.get_taxid_translator([n.name for n in t.traverse() if n.name])
    except Exception as e:
        print (e)
        sci_names = {}
        print ("NO TAXIDS")
        raise
        
    print ("Precomputing distance matrix...")
    content, node2node_dist = distance_matrix(t)

    # GRab level nodes
    target_level_nodes = set()
    node2level_name = {}
    level_sizes = []
    for lv in levels:
        print ("Preprocessing level", lv)
        try:
            if t.name == lv:
                level_node = t
            else:
                level_node = t.search_nodes(name=lv)[0]
        except:
            try:
                branch = complete_ref.search_nodes(taxid=int(lv))[0]
            except:
                print ("   Skipping taxid not found in tree", lv)
                continue
            else:
                print ("   Found level in upper branch",  branch.name)
                while len(branch.children) == 1:
                    print ("   Found one-child node",  branch.name)
                    branch = branch.children[0]                    
                try:
                    level_node = t.search_nodes(name=str(branch.taxid))[0]
                except: 
                    print ("    Skipping translated level not in target tree ", branch.name)
                    continue
                    
        print ('   Level %s is represented as %s in the target tree' %(lv, level_node.name) )        
        if len(content[level_node]) < 3:
            print ("    Skipping level with less than 3 otus", lv)
            continue
        target_level_nodes.add(level_node)
        level_node.is_level = True
        level_node.match_level = lv
        level_sizes.append(len(content[level_node]))

    LARGEST_CLADE_SIZE = float(max(level_sizes))

    
    # Process levels and detect clades        
    for node in root_node.traverse("postorder"):        
        if getattr(node, "is_level", False):
            clade_branches = set()
            print (node.name, "level")
            node.is_level=False
            for sublevel_node in node.iter_leaves(is_leaf_fn=lambda x: getattr(x, "is_level", False)):
                for c in sublevel_node.clades:
                    clade_branches.add(c)                    
            node.is_level=True
            
            clades = scan_level(node, content, node2node_dist, thr="NA", min_clade_factor=args.min_clades, dev_thr="NA", forced_clades=clade_branches, size_penalty=args.size_penalty)
            node.clades = [c[-1] for c in clades.values()]
            
            print_clade_tree(clades, node.match_level, ncbi.translate_to_names([int(node.match_level)])[0])
            
    t.write(outfile='2_bact_tree_used.nw', format=1)
