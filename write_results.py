import sys
import time

##################################################
# Functions for reporting results
##################################################

def writeTxCount(n, nAnswer, fileOut=sys.stderr): # will be used in results
    print(f'# Enrty Size {n}:', f'{nAnswer}', sep='\t', file=fileOut)
# end: def writeTxCount

def writeTxList(lstRefSeqs, fileOut): # will be used in results
    print(
        'geneName', 'name', 'chrom', 'strand', 
        'txStart', 'txEnd', 
        'cdsStart', 'cdsEnd', 
        'exonCount', 
        'exonStarts', 'exonEnds', 
        'tx_Size' , "5'UTR_Size", 'ORF_Size', "3'UTR_Size", 
        sep='\t', file=fileOut
        )
    for refSeq in lstRefSeqs:
        print(
            refSeq.sGeneSym, 
            refSeq.sRefSeqID, 
            refSeq.nChrID, 
            refSeq.sStrand, 
            refSeq.nTxStart, 
            refSeq.nTxEnd, 
            refSeq.nOrfStart, 
            refSeq.nOrfEnd, 
            refSeq.nExonCount, 
            ','.join([str(exonStart) for exonStart in refSeq.lstExonStarts]) + ',', 
            ','.join([str(exonEnd) for exonEnd in refSeq.lstExonEnds]) + ',', 
            refSeq.nExonSeqSize, 
            refSeq.nUtr5SeqSize, 
            refSeq.nOrfSeqSize, 
            refSeq.nUtr3SeqSize, 
            sep='\t', file=fileOut
            ) # end: print
# end: def writeTxList

def interimReport(lstEntrySize):
    print(time.ctime(), end=' --- ')
    writeTxCount(len(lstEntrySize), lstEntrySize[-1])
# end: def interimReport

def results(lstEntrySize, lstRefSeqs, fileOut):
    with open(fileOut, 'w') as fOut:
        for i, nAnswer in enumerate(lstEntrySize):
            writeTxCount(i + 1, nAnswer, fOut)
        # end: for i, lstRefSeqs
        print(file=fOut) # just for empty line
        writeTxList(lstRefSeqs, fOut)
# end: def results
