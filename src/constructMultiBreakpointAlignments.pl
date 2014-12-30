use strict;

## usage: perl constructMultiBreakpointAlignments.pl [infile] [minoverhang] [fileprefix] [v; optional]

my $in = $ARGV[0];
my $minoverhang = $ARGV[1];
my $prefix = $ARGV[2];
my $out = "$prefix/multibreakpointalignments.txt";
my $align = "$prefix/adjustedalignments.txt";
my $v = 0;
if (@ARGV == 4) { $v = 1; }

####

if ($v == 1) { print "VERBOSE MODE\n"; }

open(IN,$in) or die;
open(ALIGN,">$align") or die;
open(OUT,">$out") or die;
print OUT "Name\tLen\tQStart\tQEnd\tStrand\tTChr\tTStart\tTend\tQStart...\n";
my $id = -1;
my @list = ();
my @sortedlist;
my $strobenum = 0;
my %seenalignments;
my %seenalignments5;
while(my $line = <IN>) {
    chomp($line);
    my @row = split("\t",$line);
    #print "$row[5]\n";
    my $newid = $row[0];
    if($id ne $newid) {
	if ($v==1) { print "processing id = $id\n"; }
	$strobenum++;
	if($id != -1) {
	    if ($strobenum % 1000==0) {
		print "  processing read $strobenum: $list[0][0]...\n";
	    }
	    @sortedlist = sort {$a->[1] <=> $b->[1]} @list; 
	    foreach my $i (0..@sortedlist-1) {
		my @empty = ();
		&addAlignmentPaths(\@empty,\@sortedlist,$i,$v,\%seenalignments,\%seenalignments5);
	    }
	}
	@list = ();
	foreach my $key (keys %seenalignments) { delete($seenalignments{$key}); }
	$id = $newid;
	if ($v==1) { print "done processing id (newid = $id)\n";}

    }
    my @thisline = split(/\t/,$line);
    push(@list,[@thisline]);
}

## Do last longread
@sortedlist = sort {$a->[1] <=> $b->[1]} @list;
foreach my $i (0..@sortedlist-1) {
   my @empty = ();
   &addAlignmentPaths(\@empty,\@sortedlist,$i,$v,\%seenalignments,\%seenalignments5);
}

close(IN);
close(OUT);
close(ALIGN);

sub addAlignmentPaths {
    my ($cp,$sl,$index,$v,$sa,$sa5) = @_;
    my @currentpath = @{$cp};

    my @sortedlist = @{$sl};
    #my %seenalignments = %{$sa};
    push(@currentpath,$sortedlist[$index]);
    my $i = @currentpath-1;
    # if this is not the FIRST alignment, trim!
    if ($i > 0) { 
       my $leftend = $currentpath[$i-1][2];
       my $rightstart = $currentpath[$i][1];
       if($leftend > $rightstart) {
         my $overlap = $leftend-$rightstart+1;
         &adjustAlignments(\@currentpath,$i-1,$overlap,$v);
       } 
    } 
    foreach my $j ($index+1..@sortedlist-1) {
      if($currentpath[$i][2] < $sortedlist[$j][1] ||
         ($sortedlist[$j][1]-$currentpath[$i][1] >= $minoverhang &&
	       $sortedlist[$j][2]-$currentpath[$i][2] >= $minoverhang)) {
       	&addAlignmentPaths(\@currentpath,$sl,$j,$v,$sa,$sa5);
      }
    }

    print OUT "$currentpath[0][0]\t$currentpath[0][3]";
    foreach my $c (0..@currentpath-1) {
	print OUT "\t$currentpath[$c][1]\t$currentpath[$c][2]\t$currentpath[$c][4]\t$currentpath[$c][5]\t$currentpath[$c][6]\t$currentpath[$c][7]";
	
	# ONLY write alignments if it hasn't been seen before.
	# DON"T write sequence information.
	my $strtomatch = join("",@{$currentpath[$c]});	
	if (!defined($sa->{$strtomatch})) { 
	    my @noalignments = @{$currentpath[$c]};
	    @noalignments = @noalignments[0 .. 11];
	    print ALIGN join("\t",@noalignments)."\n";
	    $sa->{$strtomatch}=1;
	} else {
	    #print "SKIPPING WRITING ALIGNMENT\n";
	}
    }
    print OUT "\n";
    #if ($currentpath[0][0] eq "longread_100910"){ close(OUT); close(ALIGN); die;}
}

# % head -n 1 pacbio.machineformat.formatted.withalign
# Name    Qstart  Qend    Qlen    Strand  Chr     Tstart  Tend    #Match  #Mismatch	#Ins    #Del    Query   Ref

sub adjustAlignments {
 my ($cp,$c,$overlap,$v) = @_;
 my @currentpath = @{$cp};
 if ($v == 1) { 
     print "\nTESTING: @{$currentpath[$c]}\n@{$currentpath[$c+1]}\n";
  print "$c to $c+1: overlap $overlap\n";
  print " $c before: Query $currentpath[$c][1] $currentpath[$c][2] and Target $currentpath[$c][6] $currentpath[$c][7]\n";
  print " Alignment Info:\n\tLen=$currentpath[$c][3]\n\tMatch=$currentpath[$c][8]\n\tMismatch=$currentpath[$c][9]\n\tIns=$currentpath[$c][10]\n\tDel=$currentpath[$c][11]\n";
  print " $c+1 before: Query $currentpath[$c+1][1] $currentpath[$c+1][2] and Target $currentpath[$c+1][6] $currentpath[$c+1][7]\n";
  print " Alignment Info:\n\tLen=$currentpath[$c+1][3]\n\tMatch=$currentpath[$c+1][8]\n\tMismatch=$currentpath[$c+1][9]\n\tIns=$currentpath[$c+1][10]\n\tDel=$currentpath[$c+1][11]\n";
 }
 my $counter = 0;
 my @q = split(//,$currentpath[$c][12]);
 my @r = split(//,$currentpath[$c][13]);
 while($counter <= $overlap) {
  my $i = @q-1;
  # Assess end of alignment
  if ($q[$i] eq $r[$i]) { $currentpath[$c][8]--; } # match
  elsif ($q[$i] eq '-') { $currentpath[$c][10]--;} # INS
  elsif ($r[$i] eq '-') { $currentpath[$c][11]--;} # DEL
  else { $currentpath[$c][9]--; } # if nothing else, mismatch
  if ($q[$i] ne '-') { $counter++; $currentpath[$c][2]--; }
  if ($r[$i] ne '-') { $currentpath[$c][7]--; }
  pop(@q);
  pop(@r);
  $currentpath[$c][12] = substr($currentpath[$c][12],0,-1);
  $currentpath[$c][13] = substr($currentpath[$c][13],0,-1);
 }
 
 @q = split(//,$currentpath[$c+1][12]);
 @r = split(//,$currentpath[$c+1][13]);
 #if ($v == 1) { 
 # print "comparing for $c+1:\n".substr($currentpath[$c+1][12],0,$overlap)."\n".substr($currentpath[$c+1][13],0,$overlap)."\n";
 #}
 $counter = 0;
 while($counter <= $overlap) {
  # Assess end of alignment
  if ($q[0] eq $r[0]) { $currentpath[$c+1][8]--; } # match
  elsif ($q[0] eq '-') { $currentpath[$c+1][10]--;} # INS
  elsif ($r[0] eq '-') { $currentpath[$c+1][11]--;} # DEL
  else { $currentpath[$c+1][9]--; } # if nothing else, mismatch
  if ($q[0] ne '-') { $counter++; $currentpath[$c+1][1]++; }
  if ($r[0] ne '-') { $currentpath[$c+1][6]++; }
  shift(@q);
  shift(@r);
  $currentpath[$c+1][12] = substr($currentpath[$c+1][12],1);
  $currentpath[$c+1][13] = substr($currentpath[$c+1][13],1);
 }
 if ($v == 1) { 
  print " $c after: Query $currentpath[$c][1] $currentpath[$c][2] and Target $currentpath[$c][6] $currentpath[$c][7]\n";
  print " Alignment Info:\n\tLen=$currentpath[$c][3]\n\tMatch=$currentpath[$c][8]\n\tMismatch=$currentpath[$c][9]\n\tIns=$currentpath[$c][10]\n\tDel=$currentpath[$c][11]\n";
  print " $c+1 after: Query $currentpath[$c+1][1] $currentpath[$c+1][2] and Target $currentpath[$c+1][6] $currentpath[$c+1][7]\n";
  print " Alignment Info:\n\tLen=$currentpath[$c+1][3]\n\tMatch=$currentpath[$c+1][8]\n\tMismatch=$currentpath[$c+1][9]\n\tIns=$currentpath[$c+1][10]\n\tDel=$currentpath[$c+1][11]\n";
 print "\n\n";
 }
 if ($currentpath[$c][7]-$currentpath[$c][6] < 0 || $currentpath[$c+1][2]-$currentpath[$c+1][1] < 0) { print "ERROR TRIMMING $currentpath[$c][0]! Dying.\n"; print "Target: $currentpath[$c][7]-$currentpath[$c][6] ?<? 0\n"; print "Query:  $currentpath[$c+1][2]-$currentpath[$c+1][1] ?<? 0\n"; die; }
 if ($currentpath[$c][2] > $currentpath[$c+1][1]) { print "TRIMMING DIDN'T WORK FOR $currentpath[$c][0]! Dying.\n"; die; }
}
