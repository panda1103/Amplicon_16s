#!/usr/bin/perl
use strict;
open IN,"$ARGV[0]";
open OUT,">$ARGV[1]";

sub Sum{
	my(@arr)=@_;
	my $len=scalar @arr;
	my $total;
	foreach(0..$len){
		$total += $arr[$_];
	}
	return $total;
}

while(<IN>){
	chomp;
	my @array=split(/\t/,$_);
	my $name=shift @array;
	my $all=pop @array;
	my $s=join("\t",@array);
	if($all=~/All/){
		print OUT "$name\t$s\n";
	}else{
		my $t=Sum(@array);
		#print "$name\t@array\t$t\n";
		if($t > 0){
			print OUT "$name\t$s\n";
		}else{
			next;
		}
	}
}
