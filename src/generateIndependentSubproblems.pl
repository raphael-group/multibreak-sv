###
### generateIndependentSubproblems.pl
### Author: Anna Ritz (anna.m.ritz@gmail.com)
### Last Modified: June 2015
###
### This perl script determines independent subproblems
### (sets of clusters) that can be run in parallel through
### MultiBreak-SV.
###
### Three arguments are required:
### (1) clusters file (*.clusters)
### (2) batch input file for GASV that lists all ESPs
### (3) output directory.
###
### The output directory will be populated with clusters files.
### Each file is annotated with the subset/iteration that the
### subproblem was determined.  Pass these files with the
### --readclustersfirst option to MultiBreak-SV to run jobs
### in parallel.
###
### Use with caution - this has not been extensively tested.
use strict;

$| = 1; # force flush right away

### Three arguments are required:
### (1) clusters file (*.clusters)
### (2) batch input file for GASV that lists all ESPs
### (3) output directory.
if (@ARGV != 3) {
    print "ERROR: three arguments required.\n";
    die;
}
my $clust = $ARGV[0]; 
my $inputfile = $ARGV[1];
my $outdir = $ARGV[2];
print "ARGS ARE $clust $inputfile $outdir\n";

### Read in clusters 
open(CLUST,$clust) or die;
my %longreads;   # hash of longread -> cID
my %clusters;    # hash of cID -> list of longread
my %longreadIDs; # hash of longread -> 0 (not seen) or 1 (seen)
my %inclusters;  # hash of ESP -> 1 (in some cluster) or 0 (not in any cluster)
my %fulllines; # hash of cID -> full line
print "reading in clusters file...\n";
while(my $line = <CLUST>) {
    chomp($line);
    my @row = split(/\t/,$line);

    ## get cID; skip maximal clusters.
    my $cid = $row[0];
    if ($cid =~ m/\./) { next; }  # skip maximal clusters
    $fulllines{$cid} = $line;
    
    ## get list of ESPs; strip longreadID from ESP name.
    my @row2 = split(/,\s*/,$row[4]);
    my @longreadlist = ();
    foreach my $esp (@row2) {
	$inclusters{$esp} = 1;
	$esp =~ m/(\d+)_*/;
	my $longread = $1;
	push(@longreadlist,$longread);
	if(!defined($longreads{$longread})) { $longreads{$longread} = [()]; }
	
	# see if this is already found
	my $found = 0;
	my @array = @{$longreads{$longread}};
	foreach my $i (0..@array-1) {
	    if($array[$i] eq $cid) { $found = 1; last; }
	}

	# if not found, add it to the longread hash table.	
	if($found != 1) {
	    push(@{$longreads{$longread}},$cid);
	}

	# set longreadID to "unseen"
	$longreadIDs{$longread} = 0;
    }
    # add list of long reads to the clusters hash.
    $clusters{$cid} = [@longreadlist];
}
print "done reading clusters file.\n";
close(CLUST);

## sort the keys.
my @keys = sort {$a<=>$b} keys %longreadIDs;
shift(@keys);
my $n = @keys;
print "".(@keys)." longreads\n"; 

## initialize the seen longreads (which will be 0)
my $seen = 0;
foreach my $k (@keys) { 
    $seen+=$longreadIDs{$k}; 
}

my %thesecids; # hash of cID -> 1 if it's in this iteration.
my %shash; # 
my $num = 1;
## while we haven't seen all fragments, find the next
## indepedent subset of long reads.
while($seen < $n) {
    my @s = (); # this independent subset

    # delete last round's iteration
    foreach my $k (keys %thesecids) { delete $thesecids{$k}; }

    # get first longread
    foreach my $k (@keys) { 
	if($longreadIDs{$k} == 0) {
	    push(@s,$k);
	    $longreadIDs{$k} = 1; # seen
	    last;
	}
    }
    my $i = 0;
    ## now, find all longreads that must be considered with this longread.
    ## this includes all clusters with any longread alignment.
    while($i < @s) {
	my $longread = $s[$i];
	print "\tlongread #$i of ".(@s).": $longread\n";
	# get cids of longread
	my @cids = @{$longreads{$longread}};
	#print "\t\tcids = @cids\n";

	## for each cID that contains $longread, add the 
	## other long reads if they haven't been added already.
	foreach my $cid (@cids) {
	    $thesecids{$cid} = 1;
	    my @longreadlist = @{$clusters{$cid}};
	    #print "\t\t$cid: @longreadlist\n";
	    foreach my $longreadID (@longreadlist) {
		# check if it's already added.
		my $found = 0;
		foreach my $current (@s) {
		    if($longreadID eq $current) { $found = 1; last; }
		}

		# add longread ID if it's not in s yet.
		if($found == 0) {
		    push(@s,$longreadID);
		    $longreadIDs{$longreadID} = 1; # seen
		}
	    }
	}
	$i++;
    }
    ## By now, thesecids contains a list of all cIDs that must be considered together.
    ## Further, @s contains all longreads that must be considered together.
    ## print these statistics.
    my @cidlist =  sort {$a<=>$b} keys %thesecids;
    @s = sort{$a<=>$b} @s;
    print "SUBSET $num\n";
    print " ".(@s)." longreads connected among ".(@cidlist)." clusters\n";
    print "\tLONGREADS: @s\n";
    print "\tCLUSTERS: @cidlist\n\n";
  
    ## repopulate the hash of longreads that are present in this
    ## iteration. This will be rearranged/deleted, so this makes a 
    ## shallow copy.
    foreach my $k (keys %shash) { delete $shash{$k}; }
    foreach my $k (@s) { $shash{$k} = 1; }

    # write files. Annotate each file with this iteration.
    open(OUT,">$outdir/subset-$num.clusters") or die;
    foreach my $cid (@cidlist) {
	print OUT "$fulllines{$cid}\n";
    }
    close(OUT);

    ## increment iteration/subset.
    $num++;

    ## recalculate the number of longreads we've seen so far.
    ## the stopping criterion is when we have seen all longreads.
    $seen = 0;
    foreach my $k (@keys) { $seen+=$longreadIDs{$k}; }
}
