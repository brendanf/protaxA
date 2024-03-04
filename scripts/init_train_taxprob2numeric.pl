#!/usr/bin/perl

while (<>) {
    ($weight,$onode,$prior,$nodetype,$node,$seqid) = split;
    if ($nodetype eq "unk") {
	$node = $onode;
    }
    elsif ($nodetype =~ /^nseq1/) {
	$nodetype = (split(/,/,$nodetype))[0];
    }
    print "$nodetype $node $seqid 1 0 1.0\n";
}
