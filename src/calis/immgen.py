import warnings
import math

aas = "ACDEFGHIKLMNPQRSTVWY"
logEnrichScores = [0.127, -0.175, 0.072, 0.325, 0.380, 
                   0.110, 0.105, 0.432, -0.700, -0.036,
                   -0.570, -0.021, -0.036, -0.376, 0.168,
                   -0.537, 0.126, 0.134, 0.719, -0.012]

posImportance_9mer   = [0, 0, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0] 

pepLensSupported = [9]


def score(pepSeq=None, HLA=None):
    warnings.warn("HLA type currently ignored in the Calis Immunogenicity score: \
    Anchor regions assumed at fixed positions 2 and 9", category=PendingDeprecationWarning)
    warnings.warn("Currently, only 9-mers are supported", category=PendingDeprecationWarning)
    pepSeq = pepSeq.strip()
    if not len(pepSeq) in pepLensSupported:
        immScore = float('nan')
        return immScore
    
    immScore = 0.0
    for i, aa in enumerate(pepSeq):
        try:
            immScore += posImportance_9mer[i] * logEnrichScores[aas.index(aa)]
        except ValueError:
            raise ValueError("Amino acid symbol '%s' not supported"%(aa))
    
    return immScore