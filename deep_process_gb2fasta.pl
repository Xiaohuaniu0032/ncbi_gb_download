use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;

# data:	2021-7-20
# longfei.fu@thermofisher.com


my ($fasta,$type,$outdir);

GetOptions(
	"i:s" => \$fasta,		# fasta need to further precess          [Need]
	"t:s" => \$type,        # HBV or HCV, used to see len filter     [Need]
	"od:s" => \$outdir,     # outdir                                 [Need]
) or die "$!\n";

# first process the header
# then store all seq into its name, and use new header in the outfile

if (not defined $fasta || not defined $type || not defined $outdir){
	die "please check your args\n";
}

if (!-d $outdir){
	`mkdir -p $outdir`;
}

# -t only accept HBV or HCV
if ($type ne "HBV" and $type ne "HCV"){
	die "-t can only be [HBV|HCV]\n";
}


############################ process header ###########################
# filter seq len
my $len_cutoff;
if ($type eq "HBV"){
	$len_cutoff = 3000; # HBV ref seq NC_003977 len is 3182
}else{
	$len_cutoff = 9000; # HCV ref seq NC_004102 len is 9646
}


my @headers;
open IN, "$fasta" or die;
while (<IN>){
	chomp;
	if (/^\>/){
		my @arr = split /\|/;
		my $len = (split /\:/, $arr[3])[1];
		next if ($len < $len_cutoff); # filter len
		my $seq_name = $arr[0];
		$seq_name =~ s/^\>//;
		push @headers, $_;
	}
}
close IN;

my @headers_new;
for my $h (@headers){
	#print "old:$h\n";
	my @h = split /\|/,$h;
	my $name = $h[0];
	$name =~ s/^\>//;
	my $gt = $h[1]; # need to process
	# my $sub_gt = $h[2]; # need to process
	my $len = $h[3];
	my $country = $h[4]; # need to precess
	my $year = $h[5];

	my $new_gt = &process_gt($gt,$type);
	#next if ($new_gt eq "NA"); # skip seq with no genotype information
	my $country_new = &process_country_info($country);

	# new name format: >name|genotype:x|len:x|country:x|year:x
	my $h_new = "\>$name\|$new_gt\|$len\|$country_new\|$year";

	# skip seq with no gt
	next if $new_gt eq "genotype:NA";
	#print "new:$h_new\n";

	push @headers_new, $h_new;
}

################################ process qualified seq fasta ##########################
my $outfile = "$outdir/$type\.NCBI.fasta";
open O, ">$outfile" or die;

# first read raw fasta
my %seq;
my $seq_name;
open IN, "$fasta" or die;
while (<IN>){
	chomp;
	if (/^\>/){
		my @arr = split /\|/, $_;
		my $n = $arr[0];
		$n =~ s/^\>//;
		$seq_name = $n;
	}else{
		push @{$seq{$seq_name}}, $_;
	}
}
close IN;

#print(Dumper(\%seq));

for my $h (@headers_new){
	print O "$h\n";
	my @h = split /\|/, $h;
	my $n = $h[0];
	$n =~ s/^\>//;
	my $seq_aref = $seq{$n};
	for my $s (@{$seq_aref}){
		print O "$s\n";
	}
}
close O;


sub process_gt{
	# process genotype
	my ($gt,$type) = @_;
	my $gt_new;
	if ($type eq "HBV"){
		# HBV do not need sub genotype info
		# 7838 genotype: C
		# 7367 genotype: D
		# 1178 genotype: A1
		# 933  genotype: C1
		# 861  genotype: B/C # recombination type
		# 160  genotype C2/Ce (skip this, the right genotype is )
		# 103  subtype: Ce/C2; genotype: C
		# 122  genotype: C/D
		# 40   genotype: C/D
 		if ($gt eq "genotype:NA"){
			$gt_new = "genotype:NA";
		}elsif ($gt =~ /(genotype:[A-Z])/){
			$gt_new = $1;
		}elsif ($gt =~ /(genotype:\s[A-Z])/){
			my $a = $1;
			my @a = split /\s/, $a;
			$gt_new = "genotype:".$a[1];
		}else{
			$gt_new = "genotype:NA";
		}
	}else{
		# HCV's genotype will be 1a/2b/3c ..., no sub_gt
		# less /home/fulongfei/projects/HCVdb/combine/NCBI.HCV.fasta|grep ">"|cut -d '|' -f 2|sort|uniq -c|sort -k1nr|less
		
		# only consider these most common format
		# 26165 genotype: 1a
		# 24487 genotype: 1b
		# 10599 genotype: 3a
		# 3534  subtype: 1b
		# 2384  genotype 1b
		# 1920  genotype: 4a
		# 1713  genotype: 1b

		if ($gt =~ /(genotype:\s\d[a-z])/){
			my $a = $1;
			my @a = split /\:/, $a;
			my $gt = $a[1];
			$gt =~ s/^\s+//;
			$gt =~ s/\s+$//;
			$gt_new = "genotype:".$gt;
		}elsif ($gt =~ /(subtype:\s\d[a-z])/){
			my $a = $1;
			my @a = split /\:/, $a;
			my $gt = $a[1];
			$gt =~ s/^\s+//;
			$gt =~ s/\s+$//;
			$gt_new = "genotype:".$gt;
		}else{
			$gt_new = "genotype:NA";
		}
	}

	return $gt_new;
}

sub process_country_info{
	# process country info
	# /home/fulongfei/projects/HCVdb/combine/NCBI.HCV.fasta|grep ">"|cut -d '|' -f 5|sort|uniq -c|sort -k1nr|less
	# 62810 country:NA
	# 36454 country:USA
	# 14334 country:China
	# 13501 country:Canada
	# 11850 country:Spain:Valencia
	# 4601  country:United Kingdom

	my ($country) = @_;
	my $country_new;
	if ($country =~ /(country:[a-zA-Z]+)/){
		$country_new = $1;
	}else{
		$country_new = "NA";
	}

	return $country_new;
}