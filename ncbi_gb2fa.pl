use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;

# 2021/7/29
# longfei.fu@thermofisher.com

my ($acc_list,$virus,$outdir) = @ARGV;


if ($virus ne "HBV" and $virus ne "HCV"){
	die "virus can only be: [HBV|HCV]\n";
}

my $runsh = "$outdir/ncbi_download.sh";
open O, ">$runsh" or die;

# step 1: download genbank
my $gb = "$outdir/download.gbk";
my $gb_log = "$outdir/download.gbk.log";

print O "python3 $Bin/bin/ncbi_fasta_download.py $acc_list $gb $gb_log\n";


# step2: convert gb to fasta
my $fa = "$outdir/download.fasta";
print O "python3 $Bin/bin/parse_gb_file.v2.py -gb $gb -of $fa\n";

# step3: refine fasta file (filter len/gt NA)
print O "perl $Bin/bin/deep_process_gb2fasta.pl -i $fa -t $virus -od $outdir >$outdir/$virus\.filter.log\n";

close O;