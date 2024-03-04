#!/usr/bin/perl

$fastafile=shift;
open(FD,$fastafile);
$count=0;
while (<FD>) {
    if (/^>/) {
	($id)=split;
	$id =~ s/^>//;
	$rseqid2num{$id}=$count;
	$count++;
    }
}
close(FD);

while (<>) {
    ($weight,$onode,$prior,$nodetype,$node,$trainseq)=split;
    if ($nodetype =~ /^nseq1/) {
	($foo,$rseq)=split(/,/,$nodetype);
	if (!exists($rseqid2num{$rseq})) {
	    die ("ERROR: no fasta num for rseq '$rseq', line $_");
	}
	$nodetype = "nseq1,$rseqid2num{$rseq}";
    }
    if (!exists($rseqid2num{$trainseq})) {
	die ("ERROR: no fasta num for trainseq '$trainseq', line $_");
    }
    print "$weight $onode $prior $nodetype $node $rseqid2num{$trainseq}\n";
}

