from ete3 import Tree
import sys

if len(sys.argv) < 3:
    print ('\nMAP_DISTANCES: Transfer branch length from a source tree into a target tree\n')
    print ('  Usage: map_distances.py TARGET_TREE SOURCE_TREE\n')
    print ('  Returns: a "mapped_tree.nw" file the containing the mapped distances\n')
    sys.exit(1)
 

tgt = sys.argv[1]
src = sys.argv[2]

#tgt_tree is tree build with ncbi taxa
tgt_tree = Tree(tgt)
#srs_tree is tree build with COGs and ete3 build
src_tree = Tree(src)

print(len(tgt_tree))
print(len(src_tree))


#Find outgroup in tgt_tree, first we tried with tgt_tree.children[0] and if this node is != 1 the outgroup will be tgt_tree.children[1]
#with this outgroup we set_outgroup in the src_tree
outgroup_spcs = tgt_tree.children[0].get_leaf_names()

if len(outgroup_spcs) == 1:
    src_tree.set_outgroup(src_tree & outgroup_spcs[0].name)
else:
    out1 = tgt_tree.children[1].get_leaf_names()[0]
    src_tree.set_outgroup(src_tree & out1)
    #no entiendo porque se repite el set_outgroup¿?
    src_tree.set_outgroup(src_tree.get_common_ancestor(outgroup_spcs))
    

#diccionario nodo: set()
src_cached_content = src_tree.get_cached_content(store_attr="name", container_type=set)
#diccionario set(): nodo
content2node = dict([(frozenset(v), k) for k, v in src_cached_content.items()])


#diccionario nodo: set() para el arbol target
tgt_cached_content = tgt_tree.get_cached_content(store_attr="name", container_type=set)

for n, content in tgt_cached_content.items():
    #Pointer to parent node
    if n.up:
        #recorremos los nodos en el arbol target y el nodo parental
        # y miramos en el dictionario conten2node del src tree los nodos que contienen ese set de taxid
        target_node = content2node[frozenset(content)]
        target_up_node = content2node[frozenset(tgt_cached_content[n.up])]
        

        dist = target_node.get_distance(target_up_node)
        print (target_node.dist, dist, target_node.is_leaf(), len(content))
       
        n.dist = dist
    
tgt_tree.write(outfile="mapped_tree.nw", format=1)
