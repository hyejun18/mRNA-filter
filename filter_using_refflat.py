from collections import defaultdict
from refseq import RefSeq

##################################################
# Functions for filtering only with RefFlat
##################################################

def readRefFlat(fileIn):
    lstRefSeqs = []
    with open(fileIn, 'r') as fIn:
        for sReadLine in fIn:
            refSeq = RefSeq()
            refSeq.parse_refFlat_line(sReadLine)
            lstRefSeqs.append(refSeq)
        # end: for sReadLine
    return lstRefSeqs
# end: def readRefFlat

def filterChr(lstRefSeqs):
    # Filter the transcripts only in chr1~22, X, and Y
    lstRefSeqs = [refSeq for refSeq in lstRefSeqs if refSeq.nChrID]
    return lstRefSeqs
# end: def filterChr

def filterNM(lstRefSeqs):
    # Filter the coding transcripts (RefSeqID == NM)
    lstRefSeqs = [refSeq for refSeq in lstRefSeqs if refSeq.get_accession_prefix() == 'NM']
    return lstRefSeqs
# end: def filterNM

def filterUniqRefSeqID(lstRefSeqs):
    # Remove the transcripts that have multiple entries in the RefFlat file
    dctCountRefSeqID = defaultdict(int)
    for refSeq in lstRefSeqs:
        dctCountRefSeqID[refSeq.sRefSeqID] += 1
    # end: for refSeq
    lstUniqRefSeqID = dict(item for item in dctCountRefSeqID.items() if item[1] == 1).keys()
    lstRefSeqs = [refSeq for refSeq in lstRefSeqs if refSeq.sRefSeqID in lstUniqRefSeqID]
    return lstRefSeqs
# end: def filterUniqRefSeqID

def removeNMD(lstRefSeqs):
    # Remove the transcripts that is potential target of NMD
    lstRefSeqs = [refSeq for refSeq in lstRefSeqs if not refSeq.is_NMD_candidate()]
    return lstRefSeqs
# def removeNMD