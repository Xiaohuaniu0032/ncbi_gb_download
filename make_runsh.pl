use strict;
use warnings;
use File::Basename;

my $d = "/home/fulongfei/projects/HBVdb";
my $list_d = "/home/fulongfei/projects/HBVdb/split_list/hbv";
my @list = glob "$list_d/list_*";

for my $list (@list){
	my $bn = basename($list);
	print "python3 $d/ncbi_fasta_download.py $list $d/download/$bn\.download.xls $d/download/$bn\.log\n";
}

