from decimal import Decimal
import numpy as np

##################################################
# Global variables
##################################################

# Const.
g_MOTIF_LEN = 7
g_FOLD_CHANGE_THRESHOLD = Decimal('-0.5')
g_RELATIVE_RISK_THRESHOLD = 1
g_PSEUDO_COUNT = 0.0000001
g_PSEUDO_ARRAY = np.array([[g_PSEUDO_COUNT, 0]] * 2)

# Type of None
g_NoneType = type(None)

# Chromosomes of interest
g_lstChrs = ['chr' + str(chrNo) for chrNo in list(range(1, 23)) + ['X', 'Y']]

# Complement: A <--> T, G <--> C
# Dictionary for the Unicode code point
g_dctComplement = str.maketrans('ACGTU', 'TGCAA')

# Codons
g_sStartCodon = 'ATG'
g_tplStopCodon = ('TAA', 'TAG', 'TGA') # Using tuple because str.endswith() requires not list but tuple.
