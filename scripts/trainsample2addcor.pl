#!/usr/bin/perl

if (scalar(@ARGV) < 2) {
    die "usage: trainsample2addcor.pl classification_file seq2taxname_file\n";
}

$infile=shift;
$seq2taxfile=shift;

open(FD,$seq2taxfile);
while (<FD>) {
    ($foo1,$foo2,$seqid,$tname)=split;
    push(@seqids,$seqid);
    push(@cortnames,$tname);
}
close(FD);

open(FD,$infile);
$i=0;
while (<FD>) {
    ($nodetype,$nid,$seqid,@a)=split;
    $n2=scalar(@a);
    $n=$n2/2;
    @names=();
    @probs=();

    if ($seqid ne $seqids[$i]) {
	die "two files not in syncrony ($seqids[$i] vs $seqid), line $i+1.\n";
    }
    
    for ($k=0; $k<$n2; $k+=2) {
	push(@names,$a[$k]);
	push(@probs,$a[$k+1]);
    }
    @index = sort { $probs[$b] <=> $probs[$a] } 0..$#probs;

    if ($n < 1) {
	$tname="__";
	$prob = 0;
    }
    else {
	$tname = $names[$index[0]];
	$prob = $probs[$index[0]];
    }
    
    if ($tname eq $cortnames[$i]) {
	$cor=1;
    }
    else {
	$cor=0;
    }
    
    print "$nodetype\t$seqid\t$prob\t$cor\t$tname\t$cortnames[$i]\n";
    $i++;
}
close(FD);
