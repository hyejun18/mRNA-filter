import os
import sys
import time
from tqdm import tqdm
from global_variables import *


##################################################
# Functions for filtering with RefFlat and Fasta
##################################################

def genChrID(): # will be used in readChrSeq
    global g_lstChrs
    for i, sChrID in enumerate(g_lstChrs, 1):
        yield i, sChrID
# end: def chrID

def readChrSeq(lstRefSeqs, faDir):
    global g_lstChrs
    lstRefSeqs.sort(key=lambda x: x.nChrID)
    chrID = genChrID()
    nChrID = -1
    for refSeq in tqdm(lstRefSeqs, desc=' Parsing FASTAs'):
        # sChrSeq <- new FASTA only when refSeq.nChrID is different from currently-open chrID of FASTA
        if refSeq.nChrID != nChrID:
            nChrID, sChrID = next(chrID)
            print(time.ctime(), f'Parsing {sChrID}.fa', sep=' --- ', file=sys.stderr)
            with open(os.path.join(faDir, sChrID + '.fa'), 'r') as fIn:
                next(fIn)
                sChrSeq = fIn.read().translate({ord('\n') : None})
        # end: if refSeq.nChrID
        refSeq.parse_fasta_seq(sChrSeq)
    # end: for refSeq
    return lstRefSeqs
# end: def readChrSeq

def removeWrongOrf(lstRefSeqs):
    # Select the transcripts with
    # 1. Start codon
    # 2. Stop codon
    # 3. ORF length == Multiple of 3
    # 4. No internal stop codon
    lstRefSeqs = [
        refSeq for refSeq in lstRefSeqs if refSeq.has_start_codon()
                                       and refSeq.has_stop_codon()
                                  and (not refSeq.nOrfSeqSize % 3)
                                  and (not refSeq.is_nonsense())
    ] # end: lstRefSeqs
    return lstRefSeqs
# end: def removeWrongOrf

def filterMajorTx(lstRefSeqs):
    # Criteria for representative transcripts: the lowest RefSeqID
    lstRefSeqs.sort(key=lambda x: int(x.sRefSeqID[3:]))
    dctRefSeqs = dict()
    # key: sGeneSym
    # value: refSeq
    for refSeq in lstRefSeqs:
        if refSeq.sGeneSym not in dctRefSeqs:
            dctRefSeqs[refSeq.sGeneSym] = refSeq
        # end: if refSeq.sGeneSym
    # end: for refSeq
    return list(dctRefSeqs.values())
# end: def filterMajorTx
