#!/usr/bin/perl

# perl taxrseq2numeric.pl taxonomy refs.fa
# fasta entry must be >seqid taxname

$targetlevel = shift;
$taxfile=shift;

open(FD,$taxfile);
while (<FD>) {
    ($nid,$pid,$level,$name)=split;
    $tname2nid{$name}=$nid;
    push(@nids,$nid);
}
close(FD);

# fasta file: each entry must contain >seqid tname

$count=0;
while (<>) {
    if (/^>/) {
	($id,$tnamefull)=split;
	$tname = join(',',(split(/,/,$tnamefull))[0..($targetlevel-1)]);
        if (exists($tname2nid{$tname})) {
	    $nid=$tname2nid{$tname};
	    push(@{$rseqs{$nid}}, $count);
	}
	$count++;
    }
}

foreach $nid (@nids) {
    if (exists($rseqs{$nid})) {
	$n=scalar(@{$rseqs{$nid}});
	print "$nid\t$n\t@{$rseqs{$nid}}\n";
    }
}

