##
## Anna Ritz (aritz@cs.brown.edu, aritz@reed.edu)
## September 2015
## Post-processing script for MultiBreak-SV output.

## Takes output of MultiBreak-SV and produces a file with
## two different scores for each cluster: an SVFrequency score
## (the proportion of times that the cluster was supported by any
## multi-breakpoint-reads in the mapping over all sampled mappings)
## and an AverageProbability score (the average of all the probabilities
## of all supporting multi-breakpoint-reads).

## Note that these values do not correspond to the final probabilities
## reported in the MultiBreak-SV paper; those scores utilized the cluster
## diagrams and the alignment probabilities to compute a probability of
## each candidate cluster based on the presence AND absence of alignments.

## Further, the AverageProbability may be misleading if there is a cluster with
## many, many pairs - suppose there are 20 pairs but 3 of them are always sampled
## (avgProb=1.0) and the rest are never sampled (avgProb=0.0).  While the SVFrequency
## will be 1.0 (there are at least two reads sampled for each mapping), the average
## probability will be 3/20=0.15.

## Additionally, the location of the SV is an estimate of the breakpoint
## region as defined in the GASV file.  Thus, not all pairs of (startloc, endloc)
## values may be valid in terms of GASV clusters.

## import statements
from optparse import OptionParser
import sys

## main
def main(args):
    usage='quick-post-process.py [options]\nOutputs a file with Cluster information, SVFrequency, and AverageProbability.'
    parser = OptionParser(usage=usage)
    parser.add_option('','--clustersfile',type='string',\
                      help='GASV clusters file. Required.')
    parser.add_option('','--finalbreakpointfile',type='string',\
                      help='Final Breakpoint File. Required.')
    parser.add_option('','--finalassignmentfile',type='string',\
                      help='Final Assignments File. Required.')
    parser.add_option('','--minsupport',type='int',default=2,\
                      help='Minimum support (# of multi-breakpoint-mappings) for SVs to count towards SVFrequency.  In other words, we "count" an SV as identified in a mapping if it contains at least this many multi-breakpoint-mappings.  Default=5.')
    parser.add_option('','--outfile',type='string',default='out.txt',\
                      help='Outfile name. Default=out.txt.')
    (opts,args) = parser.parse_args()
    if opts.clustersfile == None or opts.finalbreakpointfile == None or opts.finalassignmentfile == None:
        sys.exit('Error: script requires clusters file, finalbreakpoints file, and finalassignments file.\n')

    ## read clusters file
    ## clusterinfo is indexed by cID and is a list of [svtype,startchr,startloc,endchr,endloc,pairs]
    clusterinfo = {}
    with open(opts.clustersfile) as fin:
        for line in fin:
            row = line.strip().split('\t')
            svtype = row[3]
            startchr = row[5]
            endchr = row[6]
            pairs = row[4].split(', ')

            # get startloc and endloc range
            locs = row[7].split(', ')
            starts = []
            ends = []
            for i in range(len(locs)):
                if i % 2 == 0: # start/left/x value
                    starts+= [int(locs[i])]
                else: # end/right/y value
                    ends+=[int(locs[i])]
            clusterinfo[row[0]] = [svtype,int(startchr),[min(starts),max(starts)],int(endchr),[min(ends),max(ends)],pairs]
            

    ## read finalbreakpointsfile
    ## for each cluster, count the number of mappings that have minsupport or more assignments.
    clusters2frequencies = {}
    with open(opts.finalbreakpointfile) as fin:
        for line in fin:
            if 'ClusterID' in line: # skip header
                continue
            row = line.strip().split('\t')
            numiters = int(row[1])
            ## first two columns are clusterID and numiters; thus add 2 from minsupport.
            cumulativeval = 0.0
            for i in range(opts.minsupport+2,len(row)):
                cumulativeval+=float(row[i])
            clusters2frequencies[row[0]] = cumulativeval/numiters
            

    ## read finalassignmentsfile
    ## For each assignment, we have an average assignment.
    ## Multiple assignments may contain the same discordant pairs.
    ## However, a discordant pair is selected at most ONCE for each mapping.
    ## Thus, sum all the average assignments for each discordant pair.
    pairs2avgassignments = {} # {pairs to (possibly summed) avg assignment.}
    with open(opts.finalassignmentfile) as fin:
        for line in fin:
            if 'FragmentID' in line: # skip header
                continue
            row = line.strip().split('\t')
            if row[3] == 'error': # skip error assignment
                continue
            for pair in row[3].split(','): # for each discordant pair in the assignment...
                if pair in pairs2avgassignments:
                    pairs2avgassignments[pair]+=float(row[2])
                else:
                    pairs2avgassignments[pair]=float(row[2])


    out = open(opts.outfile,'w')
    out.write('#ClusterID\tSVFreq\tAvgAssign\tSVType\tStartChr\tStartLocRange\tEndChr\tEndLocRange\tNumPairs\tPairs\n')
    for origcid in clusterinfo:
        thisinfo = clusterinfo[origcid]
        if origcid not in clusters2frequencies:
            ## this could be one of many maximal clusters for this connected component.
            ## cid refers to the connected component (no fullstops in the name).
            cid = origcid.split('.')[0]
        else:
            cid = origcid
        if cid not in clusters2frequencies:
            sys.exit('ERROR: cid="%s" (parsed to "%s") does not appear in finalbreakpoints file.' % (origcid,cid))
        svfreq = clusters2frequencies[cid]

        val = 0
        for pair in thisinfo[5]:
            if pair not in pairs2avgassignments:
                sys.exit('ERROR: end-sequence pair "%s" does not appear in finalassignments file.' % (pair))
            val+=pairs2avgassignments[pair]
        avgval = val/len(thisinfo[5])
        
        out.write('%s\t%.4f\t%.4f\t%s\t%d\t%d-%d\t%d\t%d-%d\t%d\t%s\n' % \
                  (origcid,svfreq,avgval,thisinfo[0],thisinfo[1],thisinfo[2][0],thisinfo[2][1],thisinfo[3],thisinfo[4][0],thisinfo[4][1],len(thisinfo[5]),','.join(thisinfo[5])))

    out.close()
    print 'Wrote to %s' % (opts.outfile)
    
if __name__=='__main__':
    main(sys.argv)
