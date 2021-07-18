import sys
import os
import configparser
import argparse
import glob
from Bio import SeqIO

def parse_args():
	AP = argparse.ArgumentParser("convert genbank format into fasta format")
	AP.add_argument('-gb',help='input single genbank file', dest='genbank')
	AP.add_argument('-of',help='outfile in fasta format',dest='outfile')
	
	return AP.parse_args()


def get_seq_name(seq_record):
	name = seq_record.name
	return name

def get_seq_version(seq_record):
	v = seq_record.annotations.get("sequence_version",'NA')
	return v

def get_seq_country_and_year(seq_record):
	for feature in seq_record.features:
		if feature.type == "source":
			country = feature.qualifiers.get('country','NA')
			year = feature.qualifiers.get('collection_date','NA')

	if country != 'NA':
		country_new = country[0]
	if year != 'NA':
		year_new = year[0]

	return (country_new,year_new)

def get_seq_gt(seq_record):
	'''
	gt info may in:
		1) source.organism = "HBV genotype B"
		2) source.note = "genotype: B, subgenotype: B2"

	we will check the two info
	'''
	for feature in seq_record.features:
		if feature.type == "source":
			organism = feature.qualifiers.get('organism',['NA'])
			note = feature.qualifiers.get('note',['NA'])

			#print(organism)
			#print(note)

	# check genotype
	if 'genotype' in organism[0]:
		gt = organism[0].split(' ')[-1]
	elif 'genotype' in note[0]:
		'''
		/note="genotype: B, subgenotype: B2"
		/note="genotype: C1"
		/note="genotype: D"
		/note="genotype A1: Chronic hepatitis B"
		/note="subgenotype B4"
		/note="Complete genome; genotype C; subgenotype C2."
		/note="genotype:B"
		'''
		# use reg to extract gt
		gt = note[0]
	else:
		gt = 'NA'

	# check subgenotype
	if 'subgenotype' in note[0]:
		sub_gt = note[0]
	else:
		sub_gt = 'NA'
	
	#print(gt)
	return (gt,sub_gt)


def get_seq_len(seq_record):
	seq_len = len(seq_record.seq)
	return seq_len

def get_fasta(seq_record):
	seq_arr = seq_record.format('fasta').split('\n')
	# we need to remove '>' and later use >seq_name|version|country|year|genotype|sub.gt|len to replace it
	
	'''
	need to remove first and last item

	>>> record.format('fasta').split('\n')
	['>NC_005816.1 Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', 'TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGGGGGTAA',
	'TCTGCTCTCCTGATTCAGGAGAGTTTATGGTCACTTTTGAGACAGTTATGGAAATTAAAA', '']
	'''

	del seq_arr[0]
	del seq_arr[-1]
	fa = "".join(seq_arr)

	return fa


def main():
	args = parse_args()
	infile = args.genbank
	outfile = args.outfile

	record = SeqIO.read(infile, 'genbank')
	# >seq_name|version|country|year|genotype|sub.gt|len
	seq_name = get_seq_name(record)
	seq_v = get_seq_version(record)
	seq_country,seq_year = get_seq_country_and_year(record)
	gt,sub_gt = get_seq_gt(record)
	seq_len = get_seq_len(record)
	fa = get_fasta(record)

	of = open(outfile,'w')
	#seq_h = '>%s|version:%s|%s|%s|genotype:%s|sub_gt:%s|len:%s' % (seq_name,seq_v,seq_country,seq_year,gt,sub_gt,seq_len)
	seq_h = '>%s|genotype:%s|sub_gt:%s|len:%s|country:%s|year:%s|version:%s' %(seq_name,gt,sub_gt,seq_len,seq_country,seq_year,seq_v)
	of.write(seq_h+'\n')
	of.write(fa+'\n')
	of.close()

if __name__ == "__main__":
	main()


		



