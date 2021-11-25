from global_variables import *

class RefSeq:
    def __init__(self):
        # from RefFlat
        self.__sGeneSym = 'NULL'
        self.__sRefSeqID = 'NULL'
        self.__nChrID = 0
        self.__sStrand = 'NULL'
        self.__nTxStart = 0
        self.__nTxEnd = 0
        self.__nOrfStart = 0
        self.__nOrfEnd = 0
        self.__nExonCount = 0
        self.__lstExonStarts = []
        self.__lstExonEnds = []
        self.__nExonSeqSize = 0
        self.__nOrfSeqSize = 0
        self.__nUtr5SeqSize = 0
        self.__nUtr3SeqSize = 0
        # from Fasta
        self.__sExonSeq = 'NULL'
        self.__sOrfSeq = 'NULL'
        self.__sUtr5Seq = 'NULL'
        self.__sUtr3Seq = 'NULL'
    # end: def __init__


    ##################################################
    # Getter
    ##################################################

    @property
    def sGeneSym(self):
        return self.__sGeneSym

    @property
    def sRefSeqID(self):
        return self.__sRefSeqID

    @property
    def nChrID(self):
        return self.__nChrID

    @property
    def sStrand(self):
        return self.__sStrand

    @property
    def nTxStart(self):
        return self.__nTxStart

    @property
    def nTxEnd(self):
        return self.__nTxEnd

    @property
    def nOrfStart(self):
        return self.__nOrfStart

    @property
    def nOrfEnd(self):
        return self.__nOrfEnd

    @property
    def nExonCount(self):
        return self.__nExonCount

    @property
    def lstExonStarts(self):
        return self.__lstExonStarts

    @property
    def lstExonEnds(self):
        return self.__lstExonEnds

    @property
    def nExonSeqSize(self):
        return self.__nExonSeqSize

    @property
    def nUtr5SeqSize(self):
        return self.__nUtr5SeqSize

    @property
    def nOrfSeqSize(self):
        return self.__nOrfSeqSize

    @property
    def nUtr3SeqSize(self):
        return self.__nUtr3SeqSize

    @property
    def sExonSeq(self):
        return self.__sExonSeq

    @property
    def sUtr5Seq(self):
        return self.__sUtr5Seq

    @property
    def sOrfSeq(self):
        return self.__sOrfSeq

    @property
    def sUtr3Seq(self):
        return self.__sUtr3Seq


    ##################################################
    # Setter
    ##################################################

    @sGeneSym.setter
    def sGeneSym(self, geneSym):
        self.__sGeneSym = geneSym

    @sRefSeqID.setter
    def sRefSeqID(self, refSeqID):
        # RefSeqID should be composed of two large letters, underscore, and several arabic numerals.
        assert (refSeqID[:2].isalpha()) and (refSeqID[2] == '_') and (refSeqID[3:].isdecimal()), f'invalid RefSeqID "{refSeqID}"'
        self.__sRefSeqID = refSeqID

    @nChrID.setter
    def nChrID(self, chrID):
        # Input: chrID can be both int and str such as 'chrN'
        global g_lstChrs
        if chrID.isdecimal():
            self.nChrID = int(chrID)
        elif chrID not in g_lstChrs:
            if not chrID.startswith('chr'):
                raise ValueError(f'invalid chrID "{chrID}". It should start with "chr".')
            # end: if not
            pass # init: 0
        elif chrID == 'chrX':
            self.__nChrID = 23
        elif chrID == 'chrY':
            self.__nChrID = 24
        else:
            self.__nChrID = int(chrID[3:])
    # end: def nChrID

    @sStrand.setter
    def sStrand(self, strand):
        assert (strand == '+') or (strand == '-'), f'invalid strand "{strand}"'
        self.__sStrand = strand

    @nTxStart.setter
    def nTxStart(self, txStart):
        self.__nTxStart = int(txStart)

    @nTxEnd.setter
    def nTxEnd(self, txEnd):
        self.__nTxEnd = int(txEnd)

    # orfStart, orfEnd, exonStarts, and exonEnds should be between / equal to txStart and txEnd
    def __outOfRangeError(self, position):
        raise ValueError(f'{position} of {self.sRefSeqID} in chr{self.nChrID} is out of transcript range.')

    @nOrfStart.setter
    def nOrfStart(self, orfStart):
        orfStart = int(orfStart)
        if not (self.nTxStart <= orfStart <= self.nTxEnd):
            self.__outOfRangeError(orfStart)
        # end: if not
        self.__nOrfStart = orfStart

    @nOrfEnd.setter
    def nOrfEnd(self, orfEnd):
        orfEnd = int(orfEnd)
        if not (self.nTxStart <= orfEnd <= self.nTxEnd):
            self.__outOfRangeError(orfEnd)
        # end: if not
        self.__nOrfEnd = orfEnd

    @nExonCount.setter
    def nExonCount(self, exonCount):
        self.__nExonCount = int(exonCount)

    @lstExonStarts.setter
    def lstExonStarts(self, exonStarts):
        # Input: Both comma-separated str or list are accepted.
        if isinstance(exonStarts, str):
            exonStarts = [int(exonStart) for exonStart in exonStarts.split(',') if exonStart.isdecimal()]
        # end: if isinstance
        assert isinstance(exonStarts, list)
        assert len(exonStarts) == self.nExonCount
        for exonStart in exonStarts:
            if not (self.nTxStart <= exonStart <= self.nTxEnd):
                self.__outOfRangeError(exonStart)
        # end: for exonStart
        self.__lstExonStarts = exonStarts

    @lstExonEnds.setter
    def lstExonEnds(self, exonEnds):
        # Input: Both comma-separated str or list are accepted.
        if isinstance(exonEnds, str):
            exonEnds = [int(exonEnd) for exonEnd in exonEnds.split(',') if exonEnd.isdecimal()]
        # end: if isinstance
        assert isinstance(exonEnds, list)
        assert len(exonEnds) == self.nExonCount
        for exonEnd in exonEnds:
            if not (self.nTxStart <= exonEnd <= self.nTxEnd):
                self.__outOfRangeError(exonEnd)
        # end: for exonEnd
        self.__lstExonEnds = exonEnds

    # SeqSize
    @nExonSeqSize.setter
    def nExonSeqSize(self, exonSeqSize):
        assert isinstance(exonSeqSize, int)
        assert 0 <= exonSeqSize
        assert exonSeqSize <= (self.nTxEnd - self.nTxStart)
        self.__nExonSeqSize = exonSeqSize

    @nUtr5SeqSize.setter
    def nUtr5SeqSize(self, utr5SeqSize):
        assert isinstance(utr5SeqSize, int)
        assert 0 <= utr5SeqSize 
        assert utr5SeqSize <= self.nExonSeqSize
        self.__nUtr5SeqSize = utr5SeqSize

    @nOrfSeqSize.setter
    def nOrfSeqSize(self, orfSeqSize):
        assert isinstance(orfSeqSize, int)
        assert 0 <= orfSeqSize
        assert orfSeqSize <= self.nExonSeqSize
        self.__nOrfSeqSize = orfSeqSize

    @nUtr3SeqSize.setter
    def nUtr3SeqSize(self, utr3SeqSize):
        assert isinstance(utr3SeqSize, int)
        assert 0 <= utr3SeqSize
        assert utr3SeqSize <= self.nExonSeqSize
        self.__nUtr3SeqSize = utr3SeqSize

    # Seq
    @sExonSeq.setter
    def sExonSeq(self, exonSeq):
        assert isinstance(exonSeq, str)
        assert len(exonSeq) == self.nExonSeqSize
        self.__sExonSeq = exonSeq

    @sUtr5Seq.setter
    def sUtr5Seq(self, utr5Seq):
        assert isinstance(utr5Seq, str)
        assert len(utr5Seq) == self.nUtr5SeqSize
        self.__sUtr5Seq = utr5Seq

    @sOrfSeq.setter
    def sOrfSeq(self, orfSeq):
        assert isinstance(orfSeq, str)
        assert len(orfSeq) == self.nOrfSeqSize
        self.__sOrfSeq = orfSeq

    @sUtr3Seq.setter
    def sUtr3Seq(self, utr3Seq):
        assert isinstance(utr3Seq, str)
        assert len(utr3Seq) == self.nUtr3SeqSize
        self.__sUtr3Seq = utr3Seq


    ##################################################
    # Methods for seq size
    ##################################################

    def __orfSeqSize(self): # will be used in parse_refFlat_line
        sumOfEnds = sum(nExonEnd for nExonEnd in self.lstExonEnds if self.nOrfStart < nExonEnd < self.nOrfEnd) + self.nOrfEnd
        sumOfStarts = sum(nExonStart for nExonStart in self.lstExonStarts if self.nOrfStart < nExonStart < self.nOrfEnd) + self.nOrfStart
        return sumOfEnds - sumOfStarts
    # end: def __orfSeqSize

    def __firstUtrSeqSize(self): # will be used in parse_refFlat_line
        if self.nOrfStart == self.nTxEnd: # for NR, etc.
            return self.nExonSeqSize
        sumOfEnds = sum(nExonEnd for nExonEnd in self.lstExonEnds if nExonEnd <= self.nOrfStart) + self.nOrfStart
        sumOfStarts = sum(nExonStart for nExonStart in self.lstExonStarts if nExonStart <= self.nOrfStart)
        return sumOfEnds - sumOfStarts
    # end: def __firstUtrSeqSize

    def __lastUtrSeqSize(self): # will be used in parse_refFlat_line
        if self.nOrfStart == self.nTxEnd: # for NR, etc.
            return 0
        sumOfEnds = sum(nExonEnd for nExonEnd in self.lstExonEnds if nExonEnd >= self.nOrfEnd)
        sumOfStarts = sum(nExonStart for nExonStart in self.lstExonStarts if nExonStart >= self.nOrfEnd) + self.nOrfEnd
        return sumOfEnds - sumOfStarts
    # end: def __lastUtrSeqSize

    def __lastExonSeqSize(self):
        if self.sStrand == '+':
            lastExonIndex = self.nExonCount-1
        else:
            lastExonIndex = 0
        return self.lstExonEnds[lastExonIndex] - self.lstExonStarts[lastExonIndex]
    # end: def __lastExonSeqSize


    ##################################################
    # Methods for parsing
    ##################################################

    # for RefFlat
    def parse_refFlat_line(self, sReadLine): # will be used in readRefFlat
        lstField = sReadLine.strip().split('\t')
        self.sGeneSym = lstField[0]
        self.sRefSeqID = lstField[1]
        self.nChrID = lstField[2]
        self.sStrand = lstField[3]
        self.nTxStart = lstField[4]
        self.nTxEnd = lstField[5]
        self.nOrfStart = lstField[6]
        self.nOrfEnd = lstField[7]
        self.nExonCount = lstField[8]
        self.lstExonStarts = lstField[9]
        self.lstExonEnds = lstField[10]
        self.nExonSeqSize = sum(self.lstExonEnds) - sum(self.lstExonStarts)
        self.nOrfSeqSize = self.__orfSeqSize()
        if self.sStrand == '+':
            self.nUtr5SeqSize, self.nUtr3SeqSize = self.__firstUtrSeqSize(), self.__lastUtrSeqSize()
        else: # for '-' stranded tx
            self.nUtr5SeqSize, self.nUtr3SeqSize = self.__lastUtrSeqSize(), self.__firstUtrSeqSize()
    # end: def parse_refFlat_line

    # for Fasta
    def __reverse_complement(self): # will be used in parse_fasta_seq
        self.sExonSeq = self.sExonSeq[::-1].translate(g_dctComplement)
    # end: def __reverse_complement

    def __splicing(self, sChrSeq): # will be used in parse_fasta_seq
        exonSeq = ''
        for exonStart, exonEnd in zip(self.lstExonStarts, self.lstExonEnds):
            exonSeq += sChrSeq[exonStart:exonEnd]
        # end: for exonStart, exonEnd
        return exonSeq
    # end: def __splicing

    def parse_fasta_seq(self, sChrSeq): # will be used in readChrSeq
        self.sExonSeq = self.__splicing(sChrSeq).upper()
        # end: for exonStart, exonEnd
        if self.sStrand == '-':
            self.__reverse_complement()
        # end: if self.sStrand

        nOrfStart = self.nUtr5SeqSize
        nOrfEnd = self.nExonSeqSize - self.nUtr3SeqSize

        self.sOrfSeq = self.sExonSeq[nOrfStart:nOrfEnd]
        self.sUtr5Seq = self.sExonSeq[:nOrfStart]
        self.sUtr3Seq = self.sExonSeq[nOrfEnd:]
    # end: def parse_fasta_seq


    ##################################################
    # Methods for filtering
    ##################################################

    def get_accession_prefix(self):
        sPrefix = self.sRefSeqID[:2]
        return sPrefix
    # end: def get_accession_prefix

    def has_start_codon(self):
        global g_sStartCodon
        return self.sOrfSeq.startswith(g_sStartCodon)
    # end: def has_start_codon

    def has_stop_codon(self):
        global g_tplStopCodon
        return self.sOrfSeq.endswith(g_tplStopCodon)
    # end: def has_stop_codon

    def is_nonsense(self):
        # return True if translation might be terminated in the middle of CDS
        global g_sStartCodon, g_tplStopCodon
        sOrfSeqUpper = self.sOrfSeq
        for i in range(3, len(sOrfSeqUpper) - 3, 3): # excluding first 3 nt and last 3 nt
            if sOrfSeqUpper[i:i+3] in g_tplStopCodon:
                return True
        # end: for i
        return False
    # end: def is_nonsense

    def is_NMD_candidate(self):
        return self.nUtr3SeqSize - self.__lastExonSeqSize() > 50
    # end: def is_NMD_candidate
# end: class RefSeq