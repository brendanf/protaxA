#!/usr/bin/perl

$out_scfile=shift;

# column index starting from 0, i.e. x2 and x3 refer to 3rd and 4th columns, respectively
$x2sum=0; $x2sum2;
$x3sum=0; $x3sum3;
$count=0;
while (<>) {
    push(@datrows,$_);
    ($weight, $priprob, $foo)=split;    
    @x = split(/,/,$foo);
    $nlevels=shift(@x);
    $yi=shift(@x);
    $nrows=shift(@x);
    $ncols=shift(@x);
    for ($i=0; $i<$nrows; $i++) {
	$k = $i * $ncols;
	# skip unk rows
	if ($x[$k+1]==1) {
	    $x2sum += $x[$k+2];
	    $x2sum2 += $x[$k+2] * $x[$k+2];
	    $x3sum += $x[$k+3];
	    $x3sum2 += $x[$k+3] * $x[$k+3];
	    $count++;
	}
    }    
}

$x2mean = $x2sum / $count;
$x2sd = sqrt($x2sum2/$count - $x2mean*$x2mean);
$x3mean = $x3sum / $count;
$x3sd = sqrt($x3sum2/$count - $x3mean*$x3mean);

$s = sprintf("x2mean %.6f x2sd %.6f x3mean %.6f x3sd %.6f\n",$x2mean,$x2sd,$x3mean,$x3sd);
open(OFD,">$out_scfile") or die "ERROR: cannot write to '$out_scfile'.\n$!";
print OFD "$s";
close(OFD);

foreach $dat (@datrows) {
    ($weight, $priprob, $foo)=split(/\s/,$dat);
    @x = split(/,/,$foo);
    $nlevels=shift(@x);
    $yi=shift(@x);
    $nrows=shift(@x);
    $ncols=shift(@x);
    for ($i=0; $i<$nrows; $i++) {
	$k = $i * $ncols;
	# skip unk rows
	if ($x[$k+1]==1) {
	    $x[$k+2] = sprintf("%.6f",($x[$k+2] - $x2mean)/$x2sd);
	    $x[$k+3] = sprintf("%.6f",($x[$k+3] - $x3mean)/$x3sd);
	}
    } 
    $xmat = join(',', @x);
    print "$weight\t$priprob\t$nlevels,$yi,$nrows,$ncols,$xmat\n";
}
