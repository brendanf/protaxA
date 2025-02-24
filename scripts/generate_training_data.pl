#!/usr/bin/perl

use warnings;
use strict;

my $unkname="unk";

sub readTaxonomyFile {
    # taxonomy file format: id(integer) pid(integer) level(integer) name(string) prior(float)
    
    my ($file) = @_;
    my %node2parent = ();
    my %node2children = ();
    my %node2level = ();    
    my %name2node = ();    
    my ($id, $pid, $level, $name, $prior);
    my $maxid=0;
    my %node2prior = ();
    my %unknode = ();
    my @a=();
    open(FD,$file) or die "ERROR (readTaxonomyFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($id,$pid,$level,$name,$prior) = split;	
	next if ($id eq ""); # skip empty lines
	$node2parent{$id} = $pid;
	push(@{$node2children{$pid}},$id) if ($id != $pid);
	$node2level{$id} = $level;
	$name2node{$name} = $id;
	die "ERROR (readTaxonomyFile): cannot find prior (5th column) file '$file'. $!\n" if ($prior eq "");
	$node2prior{$id} = $prior;
	@a=split(/,/,$name);
	$unknode{$id} = 1 if ($a[$#a] eq $unkname);
    } 
    close(FD);
    return (\%node2parent,\%node2children,\%node2level, \%name2node, \%node2prior, \%unknode);
}

sub readSequence2TaxonomyFile {
    my ($file, $taxname2node) = @_;
    my %seqid2taxonomy = ();
    my %node2seqids = ();

    my ($seqid,$taxname);
    open(FD, $file) or die "ERROR (readSequence2TaxonomyFile): cannot read file '$file'. $!\n";
    while (<FD>) {
	($seqid,$taxname) = split;
	next if ($seqid eq ""); # skip empty lines
	die "ERROR (readSequence2TaxonomyFile): cannot find taxonomy node for seqid: '$seqid' taxname: '$taxname'.\n" if (!exists($taxname2node->{$taxname}));
	$seqid2taxonomy{$seqid} = $taxname2node->{$taxname};
	push(@{$node2seqids{ $taxname2node->{$taxname} }}, $seqid);
    }
    close(FD);

    return (\%seqid2taxonomy,\%node2seqids);
}

sub printTaxonomy {
    my ($n2pref, $n2cref, $n2lref, $n2seqs, $n2prior) = @_;

    foreach my $i (sort {$n2lref->{$a} <=> $n2lref->{$b} || $a <=> $b} keys %$n2lref) {
	print "LEVEL $n2lref->{$i}: '$i' parent: $n2pref->{$i}, prior: $n2prior->{$i}";
	if (exists($n2cref->{$i})) {
	    print " child nodes: @{$n2cref->{$i}}";
	}
	else {
	    print " child nodes: -";
	}
	if (exists($n2seqs->{$i})) {
	    print ", seqids: @{$n2seqs->{$i}}";
	}    
	print "\n";
    }
    return (0);
}

sub readTaxonomyRseqIdFile {
    my ($filename, $rseq_separator) = @_;
    my %node2rseq = ();
    
    open(FD, $filename) or die "ERROR (readTaxonomyRseqIdFile): cannot read file '$filename'. $!\n";
    while (<FD>) {
	my ($node, $rseqs) = split;
	@{$node2rseq{$node}} = split(/$rseq_separator/, $rseqs);
    }
    close(FD);
    
    return (\%node2rseq);
}

sub addRseqs2directSeqs {
    my ($noderef,$n2rseq, $n2seqref, $s2node) = @_;
    my ($node,$rseqid);
    foreach $node (@$noderef) {
	if (exists($n2rseq->{$node})) {
	    @{$n2seqref->{$node}} = () if (!exists($n2seqref->{$node}));
	    foreach $rseqid (@{$n2rseq->{$node}}) {
		# add rseq if it's not node's direct seq (i.e. already in node's sequence list)
		push (@{$n2seqref->{$node}}, $rseqid) if (!exists($s2node->{$rseqid}) or ($s2node->{$rseqid} ne $node));
	    }
	}
    }
    return (0);
}

sub sampleIndexWithoutReplacement {
    my ($n, $maxnum) = @_;

    # returns $n non-duplicated numbers between 0...($maxnum-1)

    if ($n > $maxnum) {
	die "ERROR (sampleIndexWithoutReplacement): can't sample $n indices between 0..($maxnum-1) without replicates.\n";
    }

    my %ret = ();
    my ($i,$j,$rint,$found);
    for ($i=0; $i<$n; $i++) {
	$rint = int(rand($maxnum));
	if (!exists($ret{$rint})) {
	    $ret{$rint} = 1;
	}
	else {
	    # find next available integer, first go up and if not found then go down
	    $found = 0;
	    for ($j=$rint+1; $j<$maxnum; $j++) {
		if (!exists($ret{$j})) {
		    $ret{$j} = 1;
		    $found=1;
		    last;
		}
	    }
	    if (!$found) {
		for ($j=$rint; $j>=0; $j--) {
		    if (!exists($ret{$j})) {
			$ret{$j} = 1;
			$found=1;
			last;
		    }
		}
	    }
	}
    }
    
    return (keys %ret);
}

sub pickRandomNode {
    my ($n2cref, $n2prior) = @_;
    # pick random (probability according to taxonomy tree) child node until leaf node (== node without child nodes)
  
    my $node = 0; # start from root
    my ($i,$n,@cumprob,$r);
    my $prob = 1;
    while (exists($n2cref->{$node}) and (scalar(@{$n2cref->{$node}}) > 0)) {
	# each child to be picked has probability according to n2prior
	$n = scalar(@{$n2cref->{$node}});
	@cumprob = (0) x $n;
	for ($i=0; $i<$n; $i++) {
	    $cumprob[$i] = $n2prior->{$n2cref->{$node}->[$i]};
	    $cumprob[$i] += $cumprob[$i-1] if ($i>0);
	}	
	$r = rand($cumprob[$n-1]);
	$i=0;
	while ($cumprob[$i] < $r) {
	    $i++;
	}
	$node = $n2cref->{$node}->[$i];
	$prob = $n2prior->{$node};
    }

    return ($node, $prob);
}

sub countCandidateNodes {
    my ($node, $n2pref, $n2cref, $n2lref, $n2seqref, $n2ndref, $target_level, $min_numseqs, $min_numsisnodes) = @_;

    # count number of descendant nodes under $node which are at $target_level and have >= $min_numseqs seqs

    if (exists($n2ndref->{$node}->{$target_level}->{$min_numseqs}->{$min_numsisnodes})) {
	return ($n2ndref->{$node}->{$target_level}->{$min_numseqs}->{$min_numsisnodes});
    }
    # required information hasn't yet been calculated so let's do so and store the result for later use

    if ($n2lref->{$node} == $target_level) {
	if (exists($n2seqref->{$node}) and (scalar(@{$n2seqref->{$node}}) >= $min_numseqs) and (scalar(@{$n2cref->{$n2pref->{$node}}}) >= $min_numsisnodes)) {
	    $n2ndref->{$node}->{$target_level}->{$min_numseqs}->{$min_numsisnodes} = 1;
	    return (1);
	}
	else {
	    # don't store 0 since only existence of hash key is needed for information
	    return(0);
	}
    }

    my $n=0;
    foreach my $cnode (@{$n2cref->{$node}}) {
	$n += countCandidateNodes($cnode, $n2pref, $n2cref, $n2lref, $n2seqref, $n2ndref, $target_level, $min_numseqs, $min_numsisnodes);
    }
    $n2ndref->{$node}->{$target_level}->{$min_numseqs}->{$min_numsisnodes} = $n if ($n > 0);
    return ($n);
}

sub findNearestSameLevelNode {
    my ($node, $n2pref, $n2cref, $n2lref, $n2seqref, $n2ndref, $min_numseqs, $n2prior, $min_numsisnodes) = @_;

    # find nearest node at the same level in taxonomy with at least $min_numseqs seqs and $min_numsisnodes sister nodes
    # if many candidates, pick one randomly (probability according to taxonomy tree)

    my $target_level = $n2lref->{$node};
    my $level = $n2lref->{$node};
    my $found=0;
    my ($i,$n,$rnode,@candidates, @cumprob, $r);

    # 1) go towards root until there is a node which has a descendant at required level with >= $min_numseqs seqs
    
    while (!$found && ($level > 0)) {
	$level--;
	$node = $n2pref->{$node}; 
	$n = countCandidateNodes($node, $n2pref, $n2cref, $n2lref, $n2seqref, $n2ndref, $target_level, $min_numseqs, $min_numsisnodes);
	$found = 1 if ($n > 0);
    }

    # 2) pick descendant node from required level (randomly according to tree structure prob if many candidates)
    
    if ($found) {
	# pick one candidate at each level until target level reached
	while ($level < $target_level) {
	    $level++;

	    @candidates = ();
	    @cumprob = ();
	    $i=0;
	    foreach $rnode (@{$n2cref->{$node}}) {
		if (exists($n2ndref->{$rnode}->{$target_level}->{$min_numseqs}->{$min_numsisnodes})) {
		    push(@candidates,$rnode); 
		    push(@cumprob, $n2prior->{$rnode});
		    $cumprob[$i] += $cumprob[$i-1] if ($i>0);
		    $i++;
		}
	    }
	    $r = rand($cumprob[$i-1]);
	    $i=0;
	    while ($cumprob[$i] < $r) {
		$i++;
	    }
	    $node = $candidates[$i];
	}
    }
    
    return ($node);
}

sub trimTaxonomy {
    my ($noderef, $n2pref, $n2cref, $n2lref, $n2prior, $targetlevel) = @_;

    # remove leaf nodes whose level < targetlevel
    my ($node, %removenode, %updatepnode, $cnode, @a);
    foreach $node (@$noderef) {
	if (($n2lref->{$node} < $targetlevel) and !exists($n2cref->{$node})) {
	    $removenode{$node} = 1;
	    $updatepnode{$n2pref->{$node}}=1;

	    delete $n2pref->{$node};
	    delete $n2cref->{$node};
	    delete $n2lref->{$node};
	    delete $n2prior->{$node};
	}
    }    

    foreach $node (keys %updatepnode) {
	if (exists($n2pref->{$node})) {
	    # hasn't been removed
	    @a = @{$n2cref->{$node}};
	    @{$n2cref->{$node}} = ();
	    foreach $cnode (@a) {
		push(@{$n2cref->{$node}}, $cnode) if (!exists($removenode{$cnode}));
	    }
	}
    }

    # prior recalculation must be done starting from root and proceeding towards leaves
    my $sumprior;
    foreach $node (sort {$n2lref->{$a} <=> $n2lref->{$b}} keys %updatepnode) {
	$sumprior=0;
	if (exists($n2cref->{$node})) {
	    foreach $cnode (@{$n2cref->{$node}}) {
		$sumprior += $n2prior->{$cnode};
	    }
	    foreach $cnode (@{$n2cref->{$node}}) {
		$n2prior->{$cnode} = $n2prior->{$cnode}/$sumprior * $n2prior->{$node}; 
	    }
	}
    }
    
    # deleted nodes still exist in @$noderef but that's no problem
    
    return (0);
}

if (scalar(@ARGV) < 7) {
    die "usage: generate_training_data.pl targetlevel taxonomy_file_with_priors seqid2taxonomy_file taxonomyid2rseqid_file num_trainsamples random_seed fulloutput(yes/no) out_file (debug_mode 0/1)\n";
}

my $targetlevel = shift;
my $taxonomyfile=shift;
my $seqid2taxidfile=shift;
my $tax2ridfile=shift;
my $num_trainseqs = shift;
my $rseed = shift;
my $fulloutput = shift;
my $printdelseqs = 0;
$printdelseqs = 1 if ($fulloutput eq "yes");
my $outfile = shift;
my $debug_mode = shift;
$debug_mode = 0 if (!defined($debug_mode));

my ($n2pref,$n2cref,$n2lref,$tname2nref,$n2prior,$unknoderef) = readTaxonomyFile($taxonomyfile);
my ($s2tref, $n2sidref) = readSequence2TaxonomyFile($seqid2taxidfile,$tname2nref);
my $n2ridref = readTaxonomyRseqIdFile($tax2ridfile, ',');

# remove leaf nodes whose level < targetlevel
my @allnodes = keys %$n2pref;
trimTaxonomy(\@allnodes, $n2pref, $n2cref, $n2lref, $n2prior, $targetlevel);

printTaxonomy($n2pref, $n2cref, $n2lref, $n2sidref, $n2prior) if ($debug_mode > 1);
addRseqs2directSeqs(\@allnodes,$n2ridref, $n2sidref, $s2tref);

srand($rseed);
my @allxdat;
my %n2nd = ();
my $n2ndref = \%n2nd;

open(FD,">$outfile") or die "ERROR: cannot write to outfile '$outfile'.$!\n";

for (my $i=1; $i<=$num_trainseqs; $i++) {
    # random leaf node
    my ($node, $priprob) = pickRandomNode($n2cref, $n2prior);
    # change priprob to be conditional on parent node
    $priprob = $priprob / $n2prior->{$n2pref->{$node}};
    my %seq_dont_use = (); # sequence id flags
    my @seqs=();
    my $trainseq = "";
    my ($ri,@ri,$seqid);
    my $deleteseqs;
    print "$i) $node $priprob level $n2lref->{$node}\n" if ($debug_mode);
    # weight 1 by default for each training sample
    print FD "1 $node " . sprintf("%.6f",$priprob);

    if (exists($n2sidref->{$node}) and (scalar(@{$n2sidref->{$node}}) >= 2)) {
	# node belongs to known taxonomy with refseqs, pick one seq, rest to reference	

	@seqs = @{$n2sidref->{$node}};
	$ri = int(rand(scalar(@seqs)));	
	$trainseq = $seqs[$ri];

	if ($printdelseqs) {
	    $seq_dont_use{$trainseq} = 1;
	    $deleteseqs = join(',',sort keys %seq_dont_use);
	    print FD " nseq2 $node $trainseq $deleteseqs\n"; 
	}
	else {
	    print FD " nseq2 $node $trainseq\n";
	} 
    }
    elsif (exists($n2sidref->{$node}) and (scalar(@{$n2sidref->{$node}}) == 1)) {
	# this node represents known species with only 1 refseq
	# find nearest node with >=2 seqs: pick one seq (which can be refseq), and one refseq, other refseqs don't use

	$node = findNearestSameLevelNode($node,$n2pref,$n2cref,$n2lref,$n2sidref,$n2ndref,2,$n2prior, 0);
	# first pick one refseq (assumption: if node has at least 2 seqs, it must have at least one rseq)
	$ri = int(rand(scalar(@{$n2ridref->{$node}})));
	my $onerefseq = $n2ridref->{$node}->[$ri];

	# then pick one seq (2 in case the first is the previously picked refseq)
	@seqs = @{$n2sidref->{$node}};
	my @ri= sort {$a <=> $b} sampleIndexWithoutReplacement(2,scalar(@seqs));
	if ($seqs[$ri[0]] ne $onerefseq) {
	    $trainseq = $seqs[$ri[0]];
	}
	else {
	    $trainseq = $seqs[$ri[1]];
	}

	if ($printdelseqs) {
	    $seq_dont_use{$trainseq} = 1;	
	    foreach $seqid (@{$n2ridref->{$node}}) {
		$seq_dont_use{$seqid} = 1 if ($seqid ne $onerefseq);
	    }
	    $deleteseqs = join(',',sort keys %seq_dont_use);
	    print FD " nseq1,$onerefseq $node $trainseq $deleteseqs\n";
	}
	else {
	    print FD " nseq1,$onerefseq $node $trainseq\n";
	}
    }
    elsif (!exists($unknoderef->{$node})) {
	# this node represents known taxonomy without refseqs
	# find nearest node with >=1 seqs: pick one seq, none to reference --> all don't use
        # instead of 1, find nearest node with >=2 seqs (this guarantees along path from parent node to root it is possible to reach current node)
	$node = findNearestSameLevelNode($node,$n2pref,$n2cref,$n2lref,$n2sidref,$n2ndref,2,$n2prior, 0);
	@seqs = @{$n2sidref->{$node}};
	$ri = int(rand(scalar(@seqs)));
	$trainseq = $seqs[$ri];

	if ($printdelseqs) {
	    $seq_dont_use{$trainseq} = 1; 
	    foreach $seqid (@{$n2ridref->{$node}}) {
		$seq_dont_use{$seqid} = 1;
	    }
	    $deleteseqs = join(',',sort keys %seq_dont_use);
	    print FD " nseq0 $node $trainseq $deleteseqs\n";
	}
	else {
	    print FD " nseq0 $node $trainseq\n";
	}
    }
    else {
	# this node represents unk node
	# find nearest node with >=1 seqs: pick one seq, none to reference --> all don't use
        # instead of 1, find nearest node with >=2 seqs (this guarantees along path from parent node to root it is possible to reach current node)
	#                                      and >= 3 sisternodes (one is unk, one is mimicking unk, the rest ones are normal)
	$node = findNearestSameLevelNode($node,$n2pref,$n2cref,$n2lref,$n2sidref,$n2ndref,2,$n2prior, 3);
	@seqs = @{$n2sidref->{$node}};
	$ri = int(rand(scalar(@seqs)));
	$trainseq = $seqs[$ri];

	if ($printdelseqs) {
	    $seq_dont_use{$trainseq} = 1; 
	    foreach $seqid (@{$n2ridref->{$node}}) {
		$seq_dont_use{$seqid} = 1;
	    }
	    $deleteseqs = join(',',sort keys %seq_dont_use);
	    print FD " unk $node $trainseq $deleteseqs\n";
	}
	else {
	    print FD " unk $node $trainseq\n";
	}
    }
}
close(FD);
