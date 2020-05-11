import sys
import numpy as np
from scipy import sparse
from Bio import SeqIO
from ete3 import Tree
from collections import Counter

def seq2vector(seq):
    #return [ord(s) for s in seq]
    vector = []
    for pos in seq: 
        if pos == "-":
            vector.append(0)
        else: 
            vector.append(ord(pos.upper()))    
    return np.array(vector, dtype="int8") #sparse.lil_matrix([vector], dtype="int8")

def load_alg(fname):
    """ 
    Loads a FASTA alignment and converts it into a sparse scipy matrix. Gaps
    are converted into zeros, so they don't consume space. 

    Returns: 
      alg_matrix: alignment in sparse matrix format
      alg_index: sequence names ordered in the same way as alg_matrix rows
    """

    index = []
    alg = []
    for r in SeqIO.parse(fname, format="fasta"):
        index.append(r.id)
        alg.append(seq2vector(r.seq))
   
    named_index = {name:i for i,name in enumerate(index)}
    return np.array(alg), named_index, index


tree_file = sys.argv[1]
alg_file = sys.argv[2]
thr = float(sys.argv[3])


alg, index, i2name = load_alg(alg_file)
tree = Tree(tree_file)
tree.set_outgroup(tree.get_midpoint_outgroup())
node2content = tree.get_cached_content(store_attr="name")

for n in tree.traverse("levelorder"):  
    if n.children:         
        ch1 = n.children[0]
        ch2 = n.children[1]

        
        leaves_left = [index[name] for name in node2content[ch1]]
        leaves_right = [index[name] for name in node2content[ch2]]
        if len(leaves_left)<3 or len(leaves_right)<3:
            continue
    
        rows, cols = alg[tuple(leaves_left),:].nonzero()        
        colres_left = Counter(cols)
        cols_left = set([c for c, count in colres_left.items() if count >=2]) #>= 0.1 * len(leaves_left) ])
        
        rows, cols = alg[tuple(leaves_right),:].nonzero()        
        colres_right = Counter(cols)
        cols_right = set([c for c, count in colres_right.items() if count >=2]) #>= 0.1 * len(leaves_left) ])

        maxcols = float(max(len(colres_left), len(colres_right)))
        overlap = len(cols_left & cols_right)/maxcols
        if overlap < thr:
            #print(ch1)
            #print(ch2)
            print( overlap, ch1.dist, ch2.dist, len(leaves_left), len(leaves_right), len(cols_left), len(cols_right), 
                i2name[leaves_left[0]], i2name[leaves_right[0]],
                sep="\t") 





