import os
import sys
from Bio import Entrez
from multiprocessing import Pool
import urllib

(LIST,OUTFILE,LOGFILE) = sys.argv[1:]

# if outfile exists, first delete
if os.path.isfile(OUTFILE):
    os.remove(OUTFILE)

if os.path.isfile(LOGFILE):
    os.remove(LOGFILE)


# open outfile
of = open(OUTFILE,'w')

# open log file
log = open(LOGFILE,'w')

Entrez.email = "longfei.fu@thermofisher.com"

acce_list = []
with open(LIST,'r') as list_fh:
    for line in list_fh:
        val = line.strip()
        acce_list.append(val)


def get_fasta(acce):
    log.write("processing genome: %s\n" % acce)
    
    try:
        handle = Entrez.efetch(db="nucleotide", id=acce, rettype="gb", retmode="text")
    except urllib.error.HTTPError:
        log.write("this accession id: %s does not exist, will skip it\n" % (acce))
    else:
        fa = handle.read().strip()    
        handle.close()
        of.write(fa+'\n')
        
    '''
    if 'complete genome' in fa:
        log.write(acce + ':is a complete genome\n')
        of.write(fa+'\n')
    else:
        log.write(acce + ':is not a complete genome, will skipped\n')
    ''' 

if __name__ == "__main__":
    for i in range(len(acce_list)):
        n = acce_list[i]
        get_fasta(n)
    log.write("end processing fasta\n")



of.close()
log.close()
