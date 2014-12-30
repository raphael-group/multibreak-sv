__doc__="""Module for parsing simple readmatcher output.  Usage: python parseRm5ForDeletionsCython.py <input.m5> <output>""" 

import sys
import os
import numpy as np
import insDelHmm
import pysam

rcMap = dict(zip("ACTGNactgn-", "TGACNtgacn-"))

# transition probabilities
t_prob = np.array([[ 0.99,
                     0.000000000000000000000000000000001, 
                     0.01], 
                   [ 0.000000000000000000000000000000001, 
                     0.99,
                     0.01], 
                   [ 0.00000000005, 
                     0.00000000005, 
                     0.9999999999]]) 


# emission probabilities = ((insertion (gap), match, deletion (gap)), ...)
# the order is insertion, deletion, normal
e_prob = np.array([[.9, .05, .05], [.05, .05, .9],[.1, .8, .1]])

t_log = np.array(np.log(t_prob), dtype='f8')

e_log = np.array(np.log(e_prob), dtype='f8')

def rcSeq (seq):
    return "".join(map(lambda x: rcMap[x], seq)[::-1])


def parseRm5():
    infn , outfn = sys.argv[1], sys.argv[2]
    out = open(outfn, 'w')
    counter = 0
    currid = ""

    for line in open(infn, buffering=100000):
        values = line.rstrip("\n").split(" ")
        query_id, query_length, query_start, query_end, query_strand = values[0], int(values[1]), int(values[2]), int(values[3]), values[4]
        if currid == query_id:
            counter += 1
            continue

        currid = query_id
        target_id, target_length, target_start, target_end, target_strand = values[6], int(values[7]), int(values[8]), int(values[9]), values[10]
        alignedQuery = values[17]
        alignedTarget = values[19]
        aligned = values[18]

        score = int(values[11])
        
        if target_strand == "-":

            alignedQuery = alignedQuery[::-1]
            alignedTarget = alignedTarget[::-1]
            aligned = aligned[::-1]

        aq = alignedQuery
        at = alignedTarget
        states = insDelHmm.outputDelInsCython (aq, at, t_log, e_log)
 
        
        
        #if 1 in states or 0 in states:
        if len(states) > 1:
            insDelHmm.printSummaryFast(query_id, target_id, target_start, target_end, target_strand, states, alignedQuery, alignedTarget, out)
        counter += 1


def buildRefDict (refFile):
    refDict = {}
    seqList = []
    seqname = None
    f = open(refFile)
    for line in f.xreadlines():
        if line[0] == ">":
            print seqname
            if len(seqList) > 0:
                refDict[seqname] = "".join(seqList)
            seqname = line.strip()[1:]
            seqList = []
        else:
            seqList.append(line.strip())
    return refDict
                       
def parseSam(refDict):
    infn , outfn = sys.argv[1], sys.argv[2]
    bamfile =  pysam.Samfile( infn, "rb" )
    out = open(outfn, 'w')
    counter = 0
    currid = ""

    tidDict = {}
    for ar in bamfile.fetch():
        query_id = ar.qname
        if ar.tid in tidDict:
            target_id = tidDict[ar.tid]
        else:
            target_id  = bamfile.getrname(ar.tid)
            tidDict[ar.tid] = target_id
        if ar.is_secondary:
            counter += 1
            continue
        target_start = ar.aend-ar.alen
        target_end = ar.aend
        target_strand = "+"
        if ar.is_reverse:
            target_strand = "-"

        if currid == query_id:
            counter += 1
            continue
        currid = query_id
        alignedpositions = ar.aligned_pairs
        
        alignedQuery = buildAlignmentSequence(ar.seq, alignedpositions, 0)
        alignedTarget = buildAlignmentSequence(refDict[target_id], alignedpositions, 1)
         

        aq = alignedQuery
        at = alignedTarget
        states = insDelHmm.outputDelInsCython (aq, at, t_log, e_log)
                       
        #if 1 in states or 0 in states:
        if len(states) > 1:
            insDelHmm.printSummaryFast(query_id, target_id, target_start, target_end, target_strand, states, alignedQuery, alignedTarget, out)
        counter += 1
        

def buildAlignmentSequence (seq, aligned_pairs, index):
    charList = []
    prevposition = -1
    for aligned_pair in aligned_pairs:
        #print aligned_pair
        aln_pos = aligned_pair[index]
        if aln_pos:
            charList.append(seq[aln_pos])
        else:
            charList.append("-")
    return "".join(charList)
 
# if you need a ref file than its not a sam file
if len(sys.argv) > 3:
    refFile = sys.argv[3]
    refDict = buildRefDict(refFile)
    parseSam(refDict)
else:
    parseRm5()

