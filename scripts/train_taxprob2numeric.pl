#!/usr/bin/perl

$taxfile=shift;
open(FD,$taxfile);
while (<FD>) {
    ($nid,$pid,$level,$name)=split;
    $name2nid{$name}=$nid;
}
close(FD);

# nodetype nid seqid node1 prob1 node2 prob2 .. nodeN probN

while (<>) {
    ($nodetype,$nid,$seqid,@a)=split;
    $n2=scalar(@a);
    $n=$n2/2;
    @names=();
    @probs=();
    
    for ($i=0; $i<$n2; $i+=2) {
	if (!exists($name2nid{$a[$i]})) {
	    die "ERROR: cannot find taxnode for '$a[$i]'. Line $_";
	}
	push(@names,$name2nid{$a[$i]});
	push(@probs,$a[$i+1]);
    }
    @index = sort { $probs[$b] <=> $probs[$a] } 0..$#probs;
    
    print "$nodetype $nid $seqid $n";
    for ($i=0; $i<$n; $i++) {
	print " $names[$index[$i]] $probs[$index[$i]]";
    }
    print "\n";
}
