##
## Modified version for HMM-only, chromosome specific runs. - AR Jan 2016
## WARNING: hg38 path and minsupport (=5) are both hard-coded in.

##
## Takes PacBio machineformat (.m5) files and prepares them
## for MultiBreak-SV.  Runs GASV as a pre-processing step.
##
## Additional pre-processing steps that may be useful:
## -- Ignore ChrY alignments if female sample
## -- Remove alignments that span the centromere
## -- Remove multi-breakpoint-mappings that imply a small breakpoint (<100bp)
## -- Change the bin size for GASV: default is that every discordant pair has a separate Lmin/Lmax. The larger the bin size, the larger the region
## of uncertainty for the breakpoints (which may correspond to alignment uncertainty)
##
##
## Questions? Anna Ritz (anna.m.ritz@gmail.com)

## import statements
from optparse import OptionParser,OptionGroup
import sys
import os
import os.path
import subprocess
from bisect import bisect_left
import re
import glob

## main
def main(args):
    usage = 'M5toMBSV.py [options] <path-to-gasv> <lib-directory> <bamprefix> <chr>\n\n\tpath-to-gasv: root directory of GASV installation, \n\t\te.g., <path-to-gasv>/gasv-read-only/. \n\t\tNecessary because scripts/sortPR.bash and bin/GASV.jar\n\t\tare called in this script.\n\tlib-directory: lib/ directory from installation.\n\tbamprefix\tprefix of bam files (e.g., <pathname>/*.sorted.bam)\n\tchr\tchromosome.'
    parser = OptionParser(usage=usage)
    parser.add_option('','--prefix',default='out',type='string',\
                          help='String to prepend to all output files. Default = "out".')
    parser.add_option('','--binsize',default=None,type='int',\
  			  help='If specified, bins the the discordant pairs by Lmin/Lmax value before running GASV. The larger the bin size, the larger the region of uncertainty for the breakpoints (which may correspond to alignment uncertainty).  Increasing this binsize will merge clusters.  Optional - default behavior is not to bin the discordant pairs, and calculate Lmin/Lmax separately for each pair (using a window of +/- 300bp).')
    parser.add_option('','--experiment',default='pacbio',type='string',\
                          help='experiment label. Default = "pacbio".')
    parser.add_option('','--lambdad',default=3,type='float',\
                          help='lambda_d for MBSV (corresponds to sequence coverage). Default = 3.0.')
    parser.add_option('','--pseq',default=0.15,type='float',\
                          help='probability of an error from sequencing. Default = 0.15 (pacbio).')
    parser.add_option('','--keepconcord',action='store_true',default=False, \
			  help='Consider multi-breakpoint-mappings from reads that have some concordant alignment.  Warning: this drastically increases the number of multi-breakpoint-mappings, and is not advised for whole-genome analysis.  Default=False.')
    parser.add_option('','--outercoords',action='store_true',default=False,\
                          help='Determine discordant pairs by their start and end alignments, rather than the start and end of the variant.  These are also referred to as the outer coordinates.  Use this option when the breakpoint regions may not be precise (e.g., in highly-repetitive regions or BLASR alignments with no additional refinement).  The default behavior is to use the inner coordinates.')
    (opts, args) = parser.parse_args()
    if len(args) != 4:
        sys.exit('Error: four positional arguments required: <path-to-gasv> <lib-directory> <bamprefix> <chr>.')
    print 'Options and arguments are: ',opts,args
    gasvdir = args[0]
    libdir = args[1]
    bamprefix = args[2]
    bamchr = args[3]
    prefix = opts.prefix
    experiment = opts.experiment

    pacbiofile = '%s-empty.m5' % (prefix)
    os.system('touch %s'  % (pacbiofile))

    ## merge bamfiles with this chromosome
    mergedbam = '%s.bam' % (prefix)
    if not os.path.isfile(mergedbam):
        files = glob.glob(bamprefix)
        print 'There are %d bam files.' % (len(files))
        if len(files) > 1:
            print 'Merging bam files...'
            cmd = 'samtools merge -R chr%s %s %s' % (bamchr,mergedbam,bamprefix)
            print cmd
            os.system(cmd)
            print 'Done merging bam files.'
        else:
            print 'Extracting region from single bam file...'
            cmd = 'samtools view -b %s chr%s > %s' % (bamprefix,bamchr,mergedbam)
            print cmd
            os.system(cmd)
            print 'Done extracting region from signle bam file.'

    ## make index.
    if not os.path.isfile('%s.bai' % (mergedbam)):
        print 'Indexing merged bam file.'
        cmd = 'samtools index %s' % (mergedbam)
        print cmd
        os.system(cmd)
        print 'Done indexing bam file.'

    ## make directories
    if not os.path.isdir('%s-FormatAlignments/' % (prefix)):
	os.system('mkdir %s-FormatAlignments' % (prefix))
    if not os.path.isdir('%s-RunGASV/' % (prefix)):
	os.system('mkdir %s-RunGASV' % (prefix))
    if not os.path.isdir('%s-RunGASV/binned-esps/' % (prefix)):
	os.system('mkdir %s-RunGASV/binned-esps/' % (prefix))
    if not os.path.isdir('%s-MBSVinputs/' % (prefix)):
	os.system('mkdir %s-MBSVinputs' % (prefix))
    if not os.path.isdir('%s-MBSVinputs/cluster-subproblems/' % (prefix)):
   	os.system('mkdir %s-MBSVinputs/cluster-subproblems/' % (prefix))

    fullespfile = '%s-FormatAlignments/empty-esps.full' % (prefix)
    mapfile = '%s-FormatAlignments/mapped-names.txt' % (prefix)
    mappedname = '%s-FormatAlignments/espmapping.txt' % (prefix)

    os.system('touch %s' % (fullespfile))
    os.system('touch %s' % (mapfile))
    os.system('touch %s' % (mappedname))

    ## Run HMM if it doesn't exist
    hmmoutput = '%s-FormatAlignments/hmm-deletions.txt' % (prefix)
    if os.path.isfile(hmmoutput):
        print 'Skipping generating deletions with the HMM...file %s already exists.' % (hmmoutput)
    else:
        print 'getting HMM Regions and writing to %s' % (hmmoutput)
        reference = '/sc/orga/projects/bashia02b/GR38/hg38.fa'
        cmd = 'python %s/parseRm5ForDeletionsCython.py %s %s %s' % (libdir,mergedbam,hmmoutput,reference)                         
        print cmd                                                              
        os.system(cmd)                                                         
        print 'done getting HMM Regions and writing to %s' % (hmmoutput) 

    ## generate HMM deleted region and add them.
    espfile,mapfile,hmmdellenfile,esps = addHMMRegions(pacbiofile,fullespfile,mapfile,mappedname,prefix,libdir,opts.outercoords)
      
    ## take outer coords or inner coords
    if opts.outercoords:
        espfile = getOuterCoords(esps,prefix)
    else:
        espfile = getInnerCoords(esps,prefix)

    ## Sort ESP file
    sortedfile = sortESPfile(espfile,gasvdir)

    ## split ESPs
    gasvinfile = splitESPs(sortedfile,mapfile,'%s-RunGASV/binned-esps/' % (prefix),'%s-RunGASV/gasv.in' % (prefix),opts.binsize)

    adjustedalignmentfile = '%s-FormatAlignments/empty-adjustedalignments.txt' % (prefix)
    os.system('touch %s' % (adjustedalignmentfile))

    #sortedfile = '%s-FormatAlignments/esps-outer-coords.txt.sorted' % (prefix)
    #hmmdellenfile='%s-FormatAlignments/hmmdel-lengths.txt' % (prefix)
    #mapfile = '%s-FormatAlignments/espmapping-with-hmmdels.txt' % (prefix)

    ## make assignment file
    assignmentfile = makeAssignmentFile(prefix,sortedfile,fullespfile,adjustedalignmentfile,mapfile,hmmdellenfile,experiment)

    ## make experiment file
    experimentfile = makeExperimentFile(prefix,experiment,opts.pseq,opts.lambdad)

    gasvinfile = '%s-RunGASV/gasv.in' % (prefix)

    ## Cluster with GASV:
    clustersfile = runGASV(prefix,gasvinfile,gasvdir)
   
    ## Split clusters file into independent subproblems.
    splitClustersFile(clustersfile,'%s-MBSVinputs/cluster-subproblems/' % (prefix))

    ## Make BED file and filter the merged BAM file for alignments.
    minsupport = 5 #hard-coded for now.
    bedfile = '%s-cluster-regions.bed' % (prefix)
    clust2bed(clustersfile,minsupport,bedfile)

    filteredbam = '%s-filtered.bam' % (prefix)
    if not os.path.isfile(filteredbam):
        print 'filtering and indexing bam file from BED regions'
        cmd = 'samtools view -hb -L %s %s > %s' % (bedfile,mergedbam,filteredbam)
        print cmd
        os.system(cmd)

        cmd = 'samtools index %s' % (filteredbam)
        print cmd
        os.system(cmd)
        print 'DONE filtering and indexing bam file from BED regions.'
    print 'Filtered BAM file is %s' % (filteredbam)

    ## get reads
    readfile = '%s-cluster-reads.txt' % (prefix)
    hmmmapfile = '%s-FormatAlignments/mapped-names-from-hmm.txt' % (prefix)
    clust2reads(clustersfile,minsupport,hmmmapfile,readfile)

    print '\nDONE. Use %s, %s, and the independent subproblem cluster files in %s to run MultiBreak-SV.' % (assignmentfile,experimentfile,'%s-MBSVinputs/cluster-subproblems/' % (prefix))

    return

def addHMMRegions(m5file,fullespfile,mapfile,orignamefile,prefix,libdir,outercoords):
    print '\nADDING DELETED REGIONS FROM HMM'

    hmmoutput = '%s-FormatAlignments/hmm-deletions.txt' % (prefix)

    ## Get original names
    orignames = {}
    readnum = 0
    with open(orignamefile) as f:
        for line in f:
            if 'Name' in line: # skip header
                continue
            (val,key) = line.split()
            orignames[key] = val
            val = int(key.split('_')[1])
            readnum = max(readnum,val)
    print '%d original names' % (len(orignames))
    readnum+=1
    print 'Maxmimum read number is %d' % (readnum)
        
    ## Read current ESPs
    esps = []
    with open(fullespfile) as f:
        for line in f:
            esps += [line.split()]
    print '%d original esps.' % (len(esps))

    ## read gaps from ESPs
    gaps = {}
    with open(mapfile) as f:
        for line in f:
            row = line.strip().split()
            e = row[1].split(',')
            g = row[2].split(',')
            for i in range(len(e)):
                gaps[e[i]] = int(g[i])

    ## read HMM deletions 
    newids = 0
    oldids = 0
    newdels = []
    with open(hmmoutput) as f:
        for line in f:
            newdels += [line.split()]
    print '%d lines from hmm deletions file'%  (len(newdels))
    
    ## Convert deletions to ESPs
    ignoredchrs = set()
    newidfile = '%s-FormatAlignments/mapped-names-from-hmm.txt' % (prefix)
    out= open(newidfile,'w')
    for read_id,target_id,start_align,end_align,del_start,del_end,strand,event_type,event_size,gap_in_query,ignore1,ignore2 in newdels:
        ## Only keep deletions larger than 200bp.
        if event_type != 'D' or int(del_end)-int(del_start)<200:
            continue

        ## format chromosome
        chrom = target_id.replace('chr','')
        if chrom == 'X':
            chrom = '23'
        elif chrom == 'Y':
            chrom = '24'
        try:
            int(chrom)
        except:
            if chrom not in ignoredchrs:
                print 'WARNING: %s is not a recognized chromosome. Ignoring this ESP.' % (chrom)
                ignoredchrs.add(chrom)
            continue
            
        if read_id not in orignames:
            ## Make a new ID for this fragment.
            orignames[read_id] = 'longread_%d'  % (readnum)
            out.write('%s\t%s\n'% (read_id,orignames[read_id]))
            readnum+=1
            newids +=1
        else:
            oldids+=1

        if outercoords:
            ## Need to account for OUTER coordinates here.  This means that the breakpoints may be anywhere within the outer coordinates
            ## of the discordant pair, rather than anywhere within the inner coordinates of the discordant pair.  Thus, the breakpoints
            ## may lie "within" the aligned region, which is allowed if the alignments may span the breakpoints (especially in the case
            ## of repetitive sequences and high-error PacBio reads).
            if orignames[read_id] in gaps: # simply add querylen to gaps
                ## gap for the outer coordinates is the length of the alignment minus the event (del) plus the gap in query
                gaps[orignames[read_id]]+=(int(end_align)-int(start_align)+1)-int(event_size)+int(gap_in_query)
            else: # add esp and gap 
                esps += [[orignames[read_id],chrom,start_align,del_start,'+',chrom,del_end,end_align,'-']]
                gaps[orignames[read_id]] = (int(end_align)-int(start_align)+1)-int(event_size)+int(gap_in_query) 
        else: 
            ## NEed to account for inner coords here.  
            if orignames[read_id] in gaps: # simply add querylen to gaps                                                                                                                  
                ## gap for the outer coordinates is the length of the alignment minus the event (del) plus the gap in query                                                               
                gaps[orignames[read_id]]+=int(event_size)-int(gap_in_query)
            else: # add esp and gap                                                                                                                          
                esps += [[orignames[read_id],chrom,start_align,del_start,'+',chrom,del_end,end_align,'-']]
                gaps[orignames[read_id]] = int(event_size)-int(gap_in_query)
                
    out.close()
    print '%d esps after adding new bp calls' % (len(esps))
    print '%d new longread IDs created; %d ids seen previously' % (newids,oldids)
    print 'New mapped names from HMM calls are located in %s' % (newidfile)
   
    ## write updated files
    fullespfile = '%s-FormatAlignments/esps-with-hmmdels.full' % (prefix)
    espout = open(fullespfile,'w')
    espmapping = '%s-FormatAlignments/espmapping-with-hmmdels.txt' % (prefix)
    mapout = open(espmapping,'w')
    lengths = '%s-FormatAlignments/hmmdel-lengths.txt' % (prefix)
    lenout = open(lengths,'w')
    for esp in esps:
        espout.write('\t'.join(esp)+'\n')
        if esp[0] not in gaps:
            sys.exit('ERROR: %s does not have a gap.'% (len(esp[0])))
        mapout.write('%s\t%s\t%d\n' % (''.join(esp),esp[0],gaps[esp[0]]))
        lenout.write('%s\t%d\n'% (esp[0],int(esp[7])-int(esp[2])+gaps[esp[0]]))
    espout.close()
    mapout.close()
    lenout.close()

    print 'final files are %s and %s' % (fullespfile,espmapping)
    return fullespfile,espmapping,lengths,esps

def getOuterCoords(esps,prefix):
    print '\nCONVERTING ESPS TO OUTER COORDS'
    newespfile = '%s-FormatAlignments/esps-outer-coords.txt' % (prefix)
    out = open(newespfile,'w')
    for esp in esps:
        newline = esp[:3] + [esp[2]] + esp[4:6] + [esp[7]] + esp[7:]
        out.write('\t'.join(newline)+'\n')
    out.close()
    print 'Wrote outer coords to %s ' %(newespfile)
    return newespfile

def getInnerCoords(esps,prefix):
    print '\nCONVERTING ESPS TO INNER COORDS'
    newespfile = '%s-FormatAlignments/esps-inner-coords.txt' % (prefix)
    out = open(newespfile,'w')
    for esp in esps:
        newline = esp[:2] + [esp[3]] + esp[3:7] + [esp[6]] + esp[8:]
        out.write('\t'.join(newline)+'\n')
    out.close()
    print 'Wrote inner coords to %s' % (newespfile)
    return newespfile

def sortESPfile(espfile,gasvdir):
    ## sort ESPs
    print 'sorting ESPs'
    cmd = 'bash %s/scripts/sortPR.bash %s' % (gasvdir,espfile)
    print cmd
    os.system(cmd)
    sortedfile = '%s.sorted' % (espfile)
    print 'Sorted file is %s' % (sortedfile)
    return sortedfile

def splitESPs(espfile,mapfile,outprefix,gasvinfile,binsize):
    print '\nSPLITTING ESP FILE'
    ## read map file.
    name2gap = {}
    with open(mapfile) as fin:
        for line in fin:
            row = line.strip().split()
            names = row[1].split(',')
            gaps = row[2].split(',')
            for i in range(len(names)):
                name2gap[names[i]]=int(gaps[i])

    maxgap = max(name2gap.values())
    print 'maxgap is %d' % (maxgap)

    if binsize == None: # no bins specified; compute Lmin/Lmax for every discordant pair.
        print 'determining Lmin/Lmax for every discordant pair.'
    else:
        bins = range(0,maxgap,binsize)
        print 'Binning discordant pairs into %d bins with binsize = %d and maxgap = %d'  % (len(bins),binsize,maxgap)
        if len(bins) < 50:
            print 'bins are ',bins
            
        ## make lists of counts, filenames, and file handles for each bin.
        counts = [0 for b in bins]
        outfilenames = ['%s/intrachrom-esps-bin-%d-%d' % (outprefix,b,b+binsize-1) for b in bins]
        outs = [open(o,'w') for o in outfilenames]

    # file for bulk GASV call.
    outgasv = open(gasvinfile,'w')
    
    
    transoutfile = '%s/translocations' % (outprefix)
    transout = open(transoutfile,'w')
    transcount = 0
    ## lmin/lmax doesn't matter for translocations
    ## write line to outgasv file so translocations are accounted for.
    outgasv.write('%s\tPR\t0\t100\n' % (transoutfile))

    ## For each discordant pair, (1) put it in a bin if there are bins or (2) write a file for it and
    ## add it as a line to the outgasv file.
    with open(espfile) as fin:
        for line in fin:
            row = line.strip().split()
            name = row[0]
            gap = name2gap[name]

            ## if translocation, write and quit.
            if gap == -1:
                ## translocation
                transout.write(line)
                transcount+=1
                continue

            if binsize == None:
                # write new file that contains just this line.
                # This could be written in a simpler way (split the list),
                # but this is more clear in case bins are not None.
                outfilename = '%s/intrachrom-%s' % (outprefix,name)
                out = open(outfilename,'w')
                out.write(line)
                out.close()

                ## Write single Lmin/Lmax value to outgasv file.
                ## Note that Window of +/- 300bp is hard-coded. This is the "buffer" around the
                ## calculated gap.  TODO: make this an argument that users can change.
                outgasv.write('%s\tPR\t%d\t%d\n' % (outfilename,max(gap-300,0),gap+300))

            else:
                ## there are bins.  find the one that this discordant pair belongs in.
                ## pos is an index into outs and counts lists.
                pos = bisect_left(bins,gap)
                if pos:
                    pos-=1
                else:
                    print name,gap,pos
                    sys.exit('ERROR: %d does not fit within gaps.' % (gap))
                #print name,gap,pos
                #print outfilenames[pos]
                #sys.exit()

                # write the line to the binned file.
                outs[pos].write(line)
                counts[pos]+=1

    ## Close all files.
    transout.close()
    if binsize != None:
        # if there is a list of file handles, close them.
        for out in outs:
            out.close()
    print '%d lines written to translocations file %s'  % (transcount,transoutfile)

    ## If there are bins, write the lines to the outgasv file.
    if binsize != None:
        for i in range(len(bins)):
            if counts[i] > 0: # don't write if it's an empty file.
                print '%d lines written to binned file %s'  % (counts[i],outfilenames[i])
                outgasv.write('%s\tPR\t%d\t%d\n' % (outfilenames[i],max(bins[i]-binsize,0),bins[i]+2*binsize))
                
    outgasv.close()
    print 'wrote to file %s' % (gasvinfile)
    return gasvinfile

def makeAssignmentFile(prefix,espfile,fullespfile,alignfile,multiblockfile,newbplengths,experiment):
    print 'Arguments: ','*'.join([prefix,espfile,fullespfile,alignfile,multiblockfile,newbplengths,experiment])
    print '\nWriting Assignment File...'
    outfile = '%s-MBSVinputs/assignments.txt' % (prefix)

    ## (1) Get ESPs
    print 'getting ESPs...'
    lines = set()
    with open(espfile,'r') as infile:
        for thisrow in infile:
            row = thisrow.strip().split()
            lines.add(row[0])
    espnames = set(e.strip() for e in lines)
    multireads = {} # {longread:set(subreads)}
    hmmdelbps = set()
    for e in lines:
        matchObj = re.match(r'(longread_.*)_(.*)-(.*)',e)
        if matchObj:
            longread = matchObj.group(1)
            read1 = matchObj.group(2)
            read2 = matchObj.group(3)
        else:
            ## This is an ESP from the HMM deletion.
            hmmdelbps.add(e)
            continue
        if longread not in multireads:
            multireads[longread] = set() 
        multireads[longread].add(read1)
        multireads[longread].add(read2)
    print '  %d ESPs from %d longreads' % (len(espnames),len(multireads))
    print [a for a in espnames][0],'...'
 #   print [a for a in multireads][0],'...'

    ## (2) Get Coordinates
    print 'getting coordinates...'
    multireadcoords = {} # {longread:{subread:coords}}
    lines = []
    with open(fullespfile,'r') as infile:
        for thisrow in infile:
            row = thisrow.strip().split()
            lines += [row]
    for l in lines:
        e = l[0]
        if e not in espnames:
            continue
        matchObj = re.match(r'(.*)_(.*)-(.*)',e)
        if matchObj:
            longread = matchObj.group(1)
            read1 = matchObj.group(2)
            read2 = matchObj.group(3)
        else:
            sys.exit('Error getting coords!')
        read1coords = l[1:4]
        read2coords = l[5:8]

        if longread not in multireadcoords:
            multireadcoords[longread] = {}

        multireadcoords[longread][read1] = tuple(read1coords)
        multireadcoords[longread][read2] = tuple(read2coords)
    print '  %d (longread) tuples.' % (len(multireadcoords))
  #  print [a for a in multireadcoords][0],'...'

    ## (3) Get [err,len] of subreads
    print 'getting errors and lengths of subreads... %s' % (alignfile)
    multireadvals = {} # {longread:{subread:[err,len]}}
    lines = []
    with open(alignfile,'r') as infile:
        for thisrow in infile:
            row = thisrow.strip().split()
            lines += [row]
    lines = [[l[0]] + l[5:12] for l in lines]
    for longread,thischr,tstart,tend,match,mismatch,ins,dels in lines:
        if thischr == 'X':
            thischr = '23'
        if thischr == 'Y':
                thischr = '24'

        if longread not in multireads or longread not in multireadcoords:
            continue
        coords = (thischr,tstart,tend)

        #length = targetend-targetstart+1+ins
        length = int(tend)-int(tstart)+1+int(ins)
        # errors = mismatch + del + ins
        err = int(mismatch)+int(dels)+int(ins)

        candidatekeys = multireadcoords[longread]
        for key in candidatekeys:
            if multireadcoords[longread][key] == coords:
                if longread not in multireadvals:
                    multireadvals[longread] = {}
                multireadvals[longread][key] = [err,length]  
    print '  %d (longread) tuples with errs and lengths' % (len(multireadvals))
  #  print [(a,multireadvals[a]) for a in multireadvals][0],'...'

    ## (4) Enumerate multiblock alignments and print to outfile
    print 'enumerating multiblock alignments...'
    out = open(outfile,'w')
    out.write('LongreadID\tDiscordantPairs\tNumErrors\tBasesInAlignment(TEnd-Tstart+1+Ins)\tExperimentLabel\n')
    lines = set()
    with open(multiblockfile,'r') as infile:
        for thisrow in infile:
            row = thisrow.strip().split()
            lines.add(row[1])
    numlines = 0
    for line in lines:
        esps = line.strip().split(',')
        skip = 0
        for e in esps:
            if e not in espnames or e in hmmdelbps: ## hmm deletion; will handle later.
                skip = 1
        if skip == 1:
            continue

        ## get subreads for the multiblock alignment
        subreads = set()
        for e in esps:
            matchObj = re.match(r'(longread_.*)_(.*)-(.*)',e)
            if matchObj:
                longread = matchObj.group(1)
                read1 = matchObj.group(2)
                read2 = matchObj.group(3)
            else:
                sys.exit('Error getting subreads for %s!' % (e))
            subreads.add(read1)
            subreads.add(read2)
        if longread not in multireadvals:
            print 'ERROR! %s not in multireadvals'  %(longread)
            sys.exit()
        ## get total length and errors
        err = 0
        length = 0
        for s in subreads:
            if s not in multireadvals[longread]:
                sys.exit('ERROR! %s not in multireadvals[%s]' % (s,longread))
            err+=multireadvals[longread][s][0]
            length+=multireadvals[longread][s][1]

        ## write esps.
        espline = ','.join(esps)
        out.write('%s\t%s\t%d\t%d\t%s\n' % (longread,espline,err,length,experiment))
        numlines+=1 

    if newbplengths != '':
        lengths = {}
        with open(newbplengths,'r') as infile:
            for line in infile:
                row = line.strip().split()
                lengths[row[0]] = row[1]
        err = 0.15
        print '%d deletions from HMM breakpoints' % (len(hmmdelbps))
        for esp in hmmdelbps:
            length = int(lengths[esp])
            numerrors = int(length*0.15)
            out.write('%s\t%s\t%d\t%d\t%s\n' % (esp,esp,numerrors,length,experiment))

    print '  wrote %d lines to file %s' % (numlines,outfile)
    out.close()

    return outfile

def makeExperimentFile(prefix,experiment,pseq,lambdad):
    print '\nMaking Experiment File...'
    experimentfile = '%s-MBSVinputs/experiments.txt' % (prefix)
    outfile = open(experimentfile,'w')
    outfile.write('ExperimentLabel\tPseq\tLambdaD\n')
    outfile.write('%s\t%f\t%f\n' % (experiment,pseq,lambdad))
    outfile.close()    
    return experimentfile

def runGASV(prefix,gasvinfile,gasvdir):
    print '\nRUNNING GASV'

    # get output dir from prefix
    outputdir = '%s-RunGASV/' % (prefix)
    print 'writing to output directory',outputdir

    cmd = 'java -Xms2g -Xmx5g -jar %s/bin/GASV.jar --cluster --batch --maximal --output regions --nohead --minClusterSize 1 --outputdir %s --verbose %s'  % (gasvdir,outputdir,gasvinfile)
    print cmd
    os.system(cmd)
    clustersfile = '%s.clusters' % (gasvinfile)
    print 'Final file is %s' % (clustersfile)
    return clustersfile

def splitClustersFile(clustersfile,outdir):

    print '\nSPLITTING CLUSTERS FILE FOR MBSV...'
    frag2clust = {} # {frag: (set of clusters)}
    clust2frag = {} # {cid: (set of frags)}
    clusterlines = {} # {cid: string}

    tot = 0
    with open(clustersfile) as fin:
        for line in fin:
            tot+=1
            if tot % 10000 == 0:
                print ' line %d' % (tot)

            row = line.split('\t')
            cid = row[0]
            clusterlines[cid] = line

            esprow = row[4].split(', ')
            fragrow = set()
            ## if an esp is a longread, make it just the longreadID.
            for esp in esprow:
                # format is 'longread_2348257_0.0.1-0.0.2'
                # add just the fragment index
                line = esp.split('_')[1]
                fragrow.add(line)

 	    ## update clusters to fragment index; fragment indices to cluster.
            clust2frag[cid] = set(fragrow)
            for frag in fragrow:
                if frag not in frag2clust:
                    frag2clust[frag] = set()
                frag2clust[frag].add(cid)
        
    print '%d clusters and %d fragments' % (len(clust2frag),len(frag2clust))

    ## For each iteraction, get independent subsets.
    num=1
    while len(clust2frag) > 0:
        seenfrags = set()
        seenclusters = set()

        frag = frag2clust.keys().pop()
        seenfrags.add(frag)
        cids = frag2clust[frag]

        while len(cids.difference(seenclusters)) > 0:
            # first get new fragments from the new cids
            frags = set()
            for cid in cids.difference(seenclusters):
                frags.update(clust2frag[cid])        
            seenclusters.update(cids)

            # now get new cids from the new fragments
            newfrags = frags.difference(seenfrags)

            cids = set()
            for frag in newfrags:
                cids.update(frag2clust[frag])

            seenfrags.update(newfrags)

        # write clusters to file
        out = open(outdir+'/iter%d.clusters' % (num),'w')
        for cid in seenclusters:
            out.write(clusterlines[cid])
        out.close()

        if num % 1000 == 0:
            print ' Iteration %d' % (num)
        num+=1

        # delete entries that we've seen
        for cid in seenclusters:
            del clust2frag[cid]
        for frag in seenfrags:
            del frag2clust[frag]

    print 'Final Iteration: %d' % (num-1)
    return    

def clust2bed(clustersfile,minsupport,outfile):
    ## assumes intrachromosomal rearrangmeents only.

    out = open(outfile,'w') or die
    counter = 0
    with open(clustersfile) as fin:
        for line in fin:
            row = line.strip().split('\t')
            cid = row[0]
            support = row[1]
        ## ignore low-support clusters or non-maximal clusters.
            if '.' in cid or int(support) < minsupport:
                continue
            counter +=1
            name = '%s_support_%s' % (cid,support)
            chrom = row[6]
            coords = [int(a) for a in row[7].split(', ')]
            out.write('chr%s\t%d\t%d\t%s\n' % (chrom,min(coords),max(coords),name))
    out.close()
    print 'Wrote %d lines to %s' % (counter,outfile)
    return

def clust2reads(clustersfile,minsupport,mappingfile,outfile):
    mappedreads = {}
    with open(mappingfile) as fin:
        for line in fin:
            row = line.strip().split()
            mappedreads[row[1]] = row[0]
            
    out = open(outfile,'w') or die
    counter = 0
    with open(clustersfile) as fin:
        for line in fin:
            row = line.strip().split('\t')
            cid = row[0]
            support = row[1]
        ## ignore low-support clusters or non-maximal clusters.                                            
            if '.' in cid or int(support) < minsupport:
                continue
            counter +=1
            name = '%s_support_%s' % (cid,support)
            esps = [e for e in row[4].split(', ')]
            mappedesps = set()
            for e in esps:
                mappedesps.add(mappedreads[e])
            out.write('%s\t%s\t%s\n' % (name,support,','.join([e for e in mappedesps])))
    out.close()
    print 'Wrote %d lines to %s' % (counter,outfile)

if __name__=='__main__':
    main(sys.argv)
