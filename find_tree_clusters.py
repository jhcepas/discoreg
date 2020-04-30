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
   
    return sparse.lil_matrix(alg), index



def alg_conservation(alg, rows):
    """
    given an alignment matrix and an arbitrary set/list of rows, calculates several quality parameters: 

    Currently exploring: 
    - mean conservation per column
    - gappyness
    
    """

    con = []
    
    matrix = alg[tuple(rows),]
    size = matrix.shape[0] * matrix.shape[1]
    
    gappyness = (matrix.size / size)
        
    for col in matrix.T:
        values = col[col.nonzero()].data[0]
        counter = Counter(values)
        if counter:    
            most_common = float(counter.most_common()[0][1])
            con.append(most_common/len(values))
    con = np.array(con)
    return con.mean(), gappyness
 
# loads alg, tree
alg, alg_index = load_alg(sys.argv[2])
name2algindex = {name:i for i, name in enumerate(alg_index)}
tree = Tree(sys.argv[1])
tree.set_outgroup(tree.get_midpoint_outgroup())

print (tree)
tree.show()

# Creates a node to sequence cache 
for leaf in tree: 
    leaf.seqindex = name2algindex[leaf.name]

n2seqs = tree.get_cached_content(store_attr="seqindex")

# Iters each internal node in the tree and calculate sub-alg quality
for n in tree.traverse("level_order"):
    if n.children:
        con = alg_conservation(alg, list(n2seqs[n]))
        print(n, con)
        input()
