#!/usr/bin/perl
    
# equal prior for leaf nodes (except unks)

$unkprior_str = shift;
$unkname = "unk";

@unkprior=split(/,/,$unkprior_str);

$maxlevel=0;
%num_nodes = ();
while (<>) {
    s/\s+$//;
    ($id,$pid,$level,$name,$prior) = split; 
    next if ($id eq ""); # skip empty lines
    $dat{$id} = $_;
    $node2parent{$id} = $pid;
    push(@nids,$id);
    push(@{$node2children{$pid}},$id) if ($id != $pid);
    $node2level{$id} = $level;
    $maxlevel = $level if ($level > $maxlevel);
    @a=split(/,/,$name);
    $isunknode{$id} = 1 if ($a[$#a] eq $unkname);
} 

if ($maxlevel != scalar(@unkprior)) {
    $n = scalar(@unkprior);
    die "ERROR: taxonomy has $maxlevel levels, but you defined $n unk priors ($unkprior_str). There should be separate unk prior for each level (in one string separated by ',').\n";    
}

# calculate number of leafnodes excluding unk nodes
$totnumleafs=0;
foreach $nid (@nids) {
    if (!exists($node2children{$nid}) and !$isunknode{$nid}) {
	$pid=$nid;
	while ($pid != 0) {
	    $pid = $node2parent{$pid};
	    $num_leafnodes{$pid}++;
	}
	# for later use, set num_leafnodes = 1 for leafnode;
	$num_leafnodes{$nid} = 1;
	$totnumleafs++;
    }
}

# add leafnode count to unk nodes
foreach $nid (sort {$node2level{$b} <=> $node2level{$a}} keys %isunknode) {
    $pid = $node2parent{$nid};
    #    $num_leafnodes{$nid} = $unkprior/(1-$unkprior)*$num_leafnodes{$pid};
    $num_leafnodes{$nid} = $unkprior[$node2level{$nid}-1];
    $pid=$nid;
    while ($pid != 0) {
	$pid = $node2parent{$pid};
	$num_leafnodes{$pid} += $num_leafnodes{$nid};
    }
    $totnumleafs += $num_leafnodes{$nid};
}

#foreach $nid (sort {$a <=> $b} keys %node2level) {
#    print "$dat{$nid}\t$num_leafnodes{$nid}\n";
#}

foreach $nid (@nids) {
    $nodeprob{$nid} = $num_leafnodes{$nid}/$totnumleafs;
}

foreach $nid (@nids) {
    $f = sprintf("%.10f",$nodeprob{$nid});
    print "$dat{$nid}\t$f\n";
}
