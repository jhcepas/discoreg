import sys
import os
import argparse 
import itertools
from collections import defaultdict, OrderedDict
import numpy as np



def iter_queries(fname):
    for line in open(fname):
        fields = line.strip().split('\t')
        yield fields

# Helper functions
def motif_sum(x, y):
    min_start = min(x[0],y[0])
    max_end = max(x[1], y[1])
    return min_start, max_end

def motif_overlap(x, y):
    return max(0, min(x[1], y[1]) - max(x[0], y[0]) + 1)

def connected_components(l):
    """ it groups all connected elements in the graph, 
        returning one group for each unconnected set of connected elements """
    out = []
    while len(l)>0:
        first, *rest = l
        lf = -1
        while len(first[-1])>lf:
            lf = len(first[-1])
            rest2 = []
            for r in rest:
                if len(first[-1].intersection(r[-1]))>0:
                    first[-1] |= r[-1]
                    first[-2] |= r[-2]
                else:
                    rest2.append(r)
            rest = rest2

        out.append(first)
        l = rest
    return out


    
def main(args):
    HSP = defaultdict(list) # store High Scoring Pairs (HSPs)
    MIN_EVALUE = 0.001
    MIN_ALG_LENGTH = 10
    MIN_OVERLAP = args.min_overlap
    HITS_FILE = args.input 
    
    
    HSP_COUNTER = 0 
    # First, stores every High Score Pair (HSP) observed for each sequence, reading the all-against-all
    # BLAST matrix. Each HSP is treated as a hit in a potential MOTIF. 
    for query, group in itertools.groupby(iter_queries(HITS_FILE), lambda x: x[0]):
        for query, hit, evalue, score, length, qstart, qend, sstart, send in group:
            evalue = float(evalue)
            length, qstart, qend, sstart, send = map(int, [length, qstart, qend, sstart, send])

            if query == hit or length < MIN_ALG_LENGTH or evalue > MIN_EVALUE:
                continue
            qid = HSP_COUNTER
            HSP_COUNTER += 1 
            sid = HSP_COUNTER 
            HSP_COUNTER += 1 

            HSP[query].append([qstart, qend, qend-qstart, qid, sid])
            HSP[hit].append([sstart, send, send-sstart, sid, qid])

    # The same sequence MOTIF could produce slighltly different overlaping HSPs.
    # In the following code, we consolidate sequence specific HSPs by groups
    # those that largerly overlap. For instance, residues 20-100 from sequence A
    # could have a very good blast hit to Sequence B. Also, residues 30-95 from
    # the same Sequence A have a hit in Sequence C. Those two HSPs in sequence A
    # should be grouped by using the largest region they cover, as they are
    # expected to represent the same conserved region (MOTIF)
    all_motifs = []
    motif2seq = {}
    for seqname, hsps in HSP.items(): 
        consolidated_motifs = [] 
        # Let's see if HSPs of seqname can be merged. For this, we start picking
        # the largest HSP, and checking if any other HSP has a suffient overlap
        # to be merged with it. When merging two HSPs, we expand the sequence
        # envelope accordingly, and we update the list HSP pairs (i.e. other
        # sequences). 
        for h in sorted(hsps, reverse=True, key=lambda x: x[2]):
            merged = False
            for c in consolidated_motifs:
                length = c[1] - c[0]
                overlap = float(motif_overlap(h, c))
                if overlap / length >= MIN_OVERLAP:
                    # if the HSPs are mergable, expand the envelope of the
                    # sequence region accordingly 
                    start, end = motif_sum(h[0:2], c[0:2])
                    c[0] = start
                    c[1] = end
                    # and update the list of HSP pairs in other sequences. 
                    c[3].update([h[3], h[4]])
                    merged = True
                    break
            if not merged:
                # If not merged, we condider the HSP a new motif region in the sequence
                consolidated_motifs.append([h[0], h[1], seqname, set(h[3:5])])

        for c in consolidated_motifs:            
            all_motifs.append([set(((c[0], c[1], c[2]),)), c[3]])
       
    # At this point, "all_motifs" contain a list of motifs per sequence. Each
    # seq-motif entry provides also a set of matching regions (HSPs) in other
    # seqs, so we have a graph conecting consolidated per-sequence HSPs. The
    # connected components function returns the clusters of motifs that are
    # connected, so each cluster represents a disconnected independent group,
    # therehore a potentially independent alignable motif block in the original
    # sequences. 
    blocks  = connected_components(all_motifs)
    for o in blocks:
        print ('MOTIF BLOCK:', len(set([x[2] for x in o[0]])),
            np.mean([x[1]-x[0] for x in o[0]]),  sep="\t")

    # We can now analyze the motif/block composition of each seq, and group
    # sequences based on similar block architecture
    seq2blocks = defaultdict(set)
    for blockid, bl in enumerate(blocks):
        for blinfo in bl[0]:
            seq2blocks[blinfo[2]].add(blockid)

    # Let's analize motif architectures
    archs = [(key, tuple(sorted(values))) for key, values in seq2blocks.items()]
    for name, arch in sorted(archs, key=lambda x: x[1]):
        print(name, arch, sep="\t")

    archcounter = defaultdict(set)
    for name, arch in sorted(archs, key=lambda x: x[1]):
        archcounter[arch].add(name)
    for arch, keys in sorted(archcounter.items(), key=lambda x: len(x[1])):
        print ('Arch:', arch, len(keys), keys, sep="\t")

    # ToDo: Extract sequences from motifs

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest='input', required=True, 
                         help="All-against-all matrix file")
    parser.add_argument("--min_overlap", dest="min_overlap", default=0.75, type=float,
                        help="Min overlap between two HSPs from the same sequence to be considered the same motif and therfore be merged. ")
    args = parser.parse_args()
    main(args)
