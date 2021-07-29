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
	"i:s" => \$fasta,       # fasta need to further precess          [Need]
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
		if ($len < $len_cutoff){
			print "[Filter Reason: Length < $len_cutoff]: $_\n";
			next;
		}
		#next if ($len < $len_cutoff); # 过滤长度
		#my $seq_name = $arr[0];
		#$seq_name =~ s/^\>//;
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
	#my $gt = (split /\:/, $h[1])[1]; # 
	my $gt = $h[1]; # genotype:genotype: D
	$gt =~ s/^genotype://;
	# my $sub_gt = $h[2]; # need to process
	my $len = $h[3];
	my $country = $h[4]; # need to precess
	my $year = $h[5];

	my $new_gt = &process_gt($gt,$type);
	#print "$new_gt\n";
	if ($new_gt eq "genotype:NA"){
		my $h_raw = $h[1];
		print "[Filter Reason: Genotype is NA]: $h\n";
		next;
	}
	
	#next if ($new_gt eq "genotype:NA"); # skip

	my $country_new = &process_country_info($country);

	# new name format: >name|genotype:x|len:x|country:x|year:x
	my $h_new = "\>$name\|$new_gt\|$len\|$country_new\|$year";
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

		# for recombination type
		# less /home/fulongfei/projects/HBVdb/combine/NCBI.HBV.fasta|grep ">"|cut -d "|" -f 2|grep "/"|sort|uniq -c|sort -k1nr|less


		# HBV do not need sub genotype info
		# 7838 genotype: C
		# 7367 genotype: D
		# 1178 genotype: A1
		# 933  genotype: C1

		# for recombination type
		# 861  genotype: B/C # recombination type
		# 160  genotype C2/Ce (recombination, will skip)
		# 103  subtype: Ce/C2; genotype: C
		# 122  genotype: C/D (recombination, keep)
		# 40   genotype: C/D
 		if ($gt eq "NA"){
			$gt_new = "genotype:NA";
		}elsif ($gt =~ /genotype:\s([A-H])$/ || $gt =~ /genotype:\s([A-H])\s/){
			# genotype: C         7838
			# genotype: D         7367
			# genotype: A         5272
			$gt_new = "genotype:".$1;
			#print "$gt_new\n";
		}elsif ($gt =~ /genotype:\s([A-H])$/){
			# genotype: a         19
			$gt_new = "genotype:".uc($1);
		}elsif ($gt =~ /genotype:([A-H])$/){
			# genotype:D          80
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype:\s([A-H])\d$/){ # 限定末尾
			# genotype: A1        1178
			# genotype: C1        933
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype:([A-H])\d$/){
			# genotype:A2        9
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype\s([A-H])$/ || $gt =~/genotype\s([A-H])\;/ || $gt =~ /genotype\s([A-H])\,/ || $gt =~ /genotype\s([A-H])\s/){
			# genotype B          659
			# genotype C          510
			# genotype E          448
			# genotype D          310
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /^([A-H])$/){
			# C          390
			# B          253
			# A          84
			# F          77
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /subgenotype\s([A-H])\d$/){
			# subgenotype A1           274
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype:\s([A-H])\d[a-z]$/){
			# genotype: F1b
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype\s([A-H])\d$/){
			# genotype:genotype C1              243
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype:\s([A-H])\;/){
			# subtype: ayw2; genotype: D; PCR_primers=fwd_name: S6, rev_name: S7                 221
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype\s([A-H])\d\/([A-H])[a-z]/ || $gt =~ /genotype\s([A-H])[a-z]\/([A-H])\d/){
			# 重组型
			# genotype C2/Ce            160
			# genotype Ce/C2            XXX

			# 真实gt应该是C,参见如下:
			# subtype: C2/Ce; genotype: C         47
			# subtype: Ce/C2; genotype: C         103
			
			# 检查S1 S2是否相等
			if ($1 eq $2){
				$gt_new = "genotype:".$1;
			}else{
				$gt_new = "genotype:NA";
			}
		}elsif ($gt =~ /genotype:\s([A-H]\/[A-H])/){
			# 重组型
			# genotype: B/C         861
			# genotype: C/D         122
			$gt_new = "genotype:".$1;
		}else{
			# 其他归为NA
			$gt_new = "genotype:NA";
			#print "$gt is NA\n";
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

		if ($gt =~ /genotype:\s(\d[a-z])$/ || $gt =~ /genotype:\s(\d[A-Z])/ || $gt =~ /genotype:\s(\d[a-z])\;/){
			# genotype: 1a           26165
			# genotype: 1b           24487
			# genotype: 3a           10599
			# genotype: 1A           108
			# genotype: 1a; frequency: .00478       92
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /subtype:\s(\d[a-z])\;/){
			# subtype: 1b;
			# previous HCV antiviral treatment: naive-all; subtype: 1b; genotype: 1               3534
			# previous HCV antiviral treatment: naive-all; subtype: 1a; genotype: 1               1546
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype\s(\d[a-z])/){
			# genotype 1b            2384
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /subtype:\s([a-z])\;\sgenotype:\s(\d)/){
			# subtype: a; genotype: 1             476 seqs
			# subtype: a; genotype: 3             393 seqs
			# subtype: b; genotype: 1             208
			# subtype: b; genotype: 1
			# subtype: a; genotype: 4             133
			$gt_new = "genotype:".$2.$1;
		}elsif ($gt =~ /subtype:\s(\d[a-z])\;/){
			# subtype: 4d;
			# subtype: 4d; genotype: 4             182 seqs
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /subtype\:{0,1}\s(\d[a-z])$/){ # 限定末尾
			# can match: subtype 1b or subtype: 4d
			# genotype: subtype 1b                123
			# genotype: 4, subtype: 4d            
			$gt_new = "genotype:".$1;
		}
		elsif ($gt =~ /genotype:(\d[a-z])/){
			# genotype:1b (193)
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /genotype:\s(\d)([a-z])([a-z])/){
			# cryoglobulin positive; genotype: 2ac (192 seqs)
			$gt_new = "genotype:".$1.$2.'/'.$1.$3; # 2a/2c
		}elsif ($gt =~ /genotype:\s(\d[a-z]\/\d[a-z])$/ || $gt =~ /genotype\:{0,1}\sRF1_(\d[a-z]\/\d[a-z])$/){
			# genotype: 2k/1b         92
			# genotype: 2k/1b         26
			# genotype RF1_2k/1b      22
			$gt_new = "genotype:".$1;
		}elsif ($gt =~ /subtype\:{0,1}\s(\d[a-z]\/\d[a-z])$/){
			# subtype: 3a/2a
			$gt_new = "genotype:".$1;
		}else{
			# genotype: 1 (will skip, because it has no sub-gt info, e.g. 1a or 1b?) # 3233 seqs
			# genotype: 3 (will skip) # 1698 seqs
			# genotype: 4 (will skip) # 1090 seqs
			# genotype: 2 (will skip) # 901 seqs
			# 1 (will skip) # 899 seqs
			# 3 (will skip) # 469 seqs
			# genotype: 5 (no sub gt info)
			# genotype: 6 (299 seqs)
			# 2
			# 4 (236)
			# genotype: 6 new variant （will skip)
			# II/1b (will skip)     7
			# genotype: 2a/c (skip) 6 seqs
			# genotype: 2a/2c/2j (skip)
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
