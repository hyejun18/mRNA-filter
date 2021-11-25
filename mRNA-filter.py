import sys
import time

from arguments import *
from filter_using_refflat import *
from filter_also_using_fasta import *
from write_results import *


##################################################
# Main
##################################################

def main():
    args = getArgs()

    lstEntrySize = []

    lstRefSeqs = readRefFlat(fileIn=args.refFlat)
    lstEntrySize.append(len(lstRefSeqs))
    interimReport(lstEntrySize) # for logging

    lstRefSeqs = filterNM(filterChr(lstRefSeqs))
    lstEntrySize.append(len(lstRefSeqs))
    interimReport(lstEntrySize) # for logging

    lstRefSeqs = filterUniqRefSeqID(lstRefSeqs)
    lstEntrySize.append(len(lstRefSeqs))
    interimReport(lstEntrySize) # for logging

    lstRefSeqs = removeNMD(lstRefSeqs)
    lstEntrySize.append(len(lstRefSeqs))
    interimReport(lstEntrySize)

    lstRefSeqs = readChrSeq(lstRefSeqs, faDir=args.fa_dir)
    lstRefSeqs = removeWrongOrf(lstRefSeqs)
    lstEntrySize.append(len(lstRefSeqs))
    interimReport(lstEntrySize) # for logging

    lstRefSeqs = filterMajorTx(lstRefSeqs)
    lstEntrySize.append(len(lstRefSeqs))
    interimReport(lstEntrySize) # for logging

    results(lstEntrySize, lstRefSeqs, fileOut=args.out_file)
# end: def main

if __name__ == '__main__':
    print(time.ctime(), 'Start', sep=' --- ', file=sys.stderr)
    main()
    print(time.ctime(), 'End', sep=' --- ', file=sys.stderr)
