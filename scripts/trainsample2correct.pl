#!/usr/bin/perl

if (scalar(@ARGV) < 2) {
    die "usage: trainsample2correct.pl level taxonomy_file trainsamplefile\n";
}

$targetlevel=shift;
$taxonomyfile=shift;

$unkname = "unk";

open(FD,$taxonomyfile);
while (<FD>) {
    ($id,$pid,$level,$name)=split;
    $node2tname{$id}=$name;
    $node2parent{$id} = $pid;
    $node2level{$id} = $level;
    $tname2node{$name}=$id;
}
close(FD);

while (<>) {
    ($weight,$onode,$priprob,$nodeinfo,$rnode,$trainseq)=split;
    if (!exists($node2tname{$rnode})) {
	die "ERROR: cannot find taxonomy node for $_";
    }    

    if ($node2level{$rnode} == $targetlevel) {
	if (($nodeinfo eq "nseq0") or ($nodeinfo eq "unk")) {
	    if ($targetlevel == 1) {
		$cor = $unkname;
	    }
	    else {
		$cor = $node2tname{$node2parent{$rnode}} . ",$unkname";
	    }
	}
	else {
	    $cor = $node2tname{$rnode};
	}
    }
    else {
	$nodeinfo = "nseq2";
	while ($node2level{$rnode} > $targetlevel) {
	    $rnode = $node2parent{$rnode};	    
	}
	$cor = $node2tname{$rnode};
    }
    
    print "$nodeinfo $rnode $trainseq $cor\n";
}
