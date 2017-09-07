import inspect
import os
import subprocess
import shutil
import io_evd.zipbin
import protbin.cls
import numpy
import warnings
import datetime

class Struct:
    def __init__(self):
        pass 


class Proteome:

    rawFields = ("pos", "AA", "C", "score", "Ident")
    dataSubDir = "data"
    protFastaSubDir = "protFasta"
    
    
    def __init__(self, netchopEnv="/tbb/local/ii/src/netchop-3.1", binName="/tbb/local/ii/src/netchop-3.1/bin/netChop", version="Cterm3_0", dataDir=None, protName=None):
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Proteome.dataSubDir))
        
        self._netchopEnv = netchopEnv
        self._binName = binName
        self._version  = version
        self._dataDir = dataDir
        self._setup_prot(protName)
    
    def _setup_prot(self, protName=None):
        if protName == None:
            raise TypeError("'protName' parameter must be specified")
        self._protName = protName
        self._protFastaDir = os.path.join(self._dataDir, Proteome.protFastaSubDir)
        self._protFastaFile = os.path.join(self._protFastaDir, "%s.fasta"%(protName))
        self._chopDataDir = os.path.join(self._dataDir, self._version)
        self._chopDataFile = os.path.join(self._chopDataDir, "%s.dat"%(self._protName))
        if not os.path.isdir(self._chopDataDir):
            os.mkdir(self._chopDataDir)
        
        if not (os.path.isfile(self._chopDataFile) or os.path.isfile(self._protFastaFile)):
            raise IOError("Proteome '%s' not recognized. Please ensure that the following file exists:\n%s"%(self._protName, self._protFastaFile))
        
    def get_scores(self):
        if not os.path.isfile(self._chopDataFile):
            dataDict = self.create_data_dict()
        else:
            inStream = io_evd.zipbin.Stream(self._chopDataFile)
            dataDict = inStream.read()
        return dataDict
    
        
    def create_data_dict(self):
        if not os.path.isfile(self._protFastaFile):
            raise IOError("Proteome '%s' not recognized. Please ensure that the following file exists:\n%s"%(self._protName, self._protFastaFile))
        
        curr_time = datetime.datetime.now();
        curr_time_str = curr_time.strftime("%Y_%m_%d_%H_%M_%S_%f")
        tempDir = os.path.join(self._chopDataDir, "%s_temp_%s"%(self._protName, curr_time_str))
        tempChopTxt = os.path.join(self._chopDataDir, "%s_temp_%s.txt"%(self._protName, curr_time_str))
        
        if not os.path.isdir(tempDir):
            os.mkdir(tempDir)    
        
        if not os.path.isfile(tempChopTxt):
            args = [self._binName, '-tdir', tempDir, self._protFastaFile]
            my_env = os.environ.copy()
            my_env["NETCHOP"] = self._netchopEnv
            with open(tempChopTxt, "wb") as f:
                proc = subprocess.Popen(args, env=my_env, stdout=f)
                proc.wait()

        protIDs = self._get_prot_keys_from_fasta()
        dataDict = netchop_raw_file_2_Dict(tempChopTxt, protIDs)
        
        outStream = io_evd.zipbin.Stream(self._chopDataFile)
        outStream.write(dataDict)
        
        shutil.rmtree(tempDir, ignore_errors=False)
        os.remove(tempChopTxt)

        return dataDict
        
    
    def _get_prot_keys_from_fasta(self):
        protIDs = []
        
        with open(self._protFastaFile, 'rb') as f:
            for line in f:
                if line[0] == ">":
                    protIDs.append(line[1:].strip())
#                     protSeqIDs.append(line[1:line.find("|",4)+1].strip()) 
#                     geneStart = line.find("GN=")+3
#                     geneEnd = line.find("PE=")
#                     geneIDs.append(line[geneStart:geneEnd].strip())
        return protIDs
    
    def _get_prot_keys_from_fasta_old(self):
        protSeqIDs = []
        geneIDs = []
        
        with open(self._protFastaFile, 'rb') as f:
            for line in f:
                if line[0] == ">":
                    protSeqIDs.append(line[1:line.find("|",4)+1].strip()) 
                    geneStart = line.find("GN=")+3
                    geneEnd = line.find("PE=")
                    geneIDs.append(line[geneStart:geneEnd].strip())
        return protSeqIDs, geneIDs   
    
    def convert_oldDat_2_newDat(self, outFile):  
        if not os.path.isfile(self._chopDataFile):
            raise IOError("Protein dat file does not exist")
        
        outFile = os.path.join(self._chopDataDir, outFile)
        
        inStream = io_evd.zipbin.Stream(self._chopDataFile)
        dataDict = inStream.read()
        
        newProtIDs = self._get_prot_keys_from_fasta()
        
        newDict = dict()
        
        for key in dataDict.keys():
            newKeyFound = False
            for newKey in newProtIDs:
                if key[0] in newKey and key[1] in newKey:
                    newDict[newKey] = dataDict[key]
                    newKeyFound = True
                    break
            if not newKeyFound:
                print "New key not found for (%s, %s)" % (key[0], key[1])

        outStream = io_evd.zipbin.Stream(outFile)
        outStream.write(newDict)
        
        
class Peptide:
    dataSubDir = "data"
    protFastaSubDir = "protFasta"
    tieBreaks       = ["MAX", "MIN", "MEAN"]
    aaExtendLeft    = 8
    aaExtendRight   = 8
    
    def __init__(self, \
                 netchopEnv = "/tbb/local/ii/src/netchop-3.1", \
                 binName    = "/tbb/local/ii/src/netchop-3.1/bin/netChop", \
                 version    = "Cterm3_0", \
                 dataDir    = None, \
                 proteomeName   = None, \
                 tieBreak   = "MAX", \
                 ):      
          
        
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Peptide.dataSubDir))
        
        self._netchopEnv = netchopEnv
        self._binName = binName
        self._version  = version
        self._dataDir = dataDir
        self._tieBreak = tieBreak
        self._setup_proteome(proteomeName)
        
        
    def _setup_proteome(self, proteomeName=None):
        if proteomeName == None:
            raise TypeError("'proteomeName' parameter not specified")
        self._proteomeName = proteomeName
        self._proteomeLoaded = False
        self._protFastaDir = os.path.join(self._dataDir, Peptide.protFastaSubDir)
        self._protFastaFile = os.path.join(self._protFastaDir, "%s.fasta"%(proteomeName))
        curr_time = datetime.datetime.now();
        curr_time_str = curr_time.strftime("%Y_%m_%d_%H_%M_%S_%f")
        self._netChopTemp = os.path.join(self._protFastaDir, "%s_netchop_evd_temp_%s"%(proteomeName, curr_time_str))
        
    def _load_proteome(self):
        if not self._proteomeLoaded:
            protBinStream = protbin.cls.Proteome(dataDir=self._dataDir, protName=self._proteomeName)
            self._proteomeDict = protBinStream.get_prots()
            self._proteomeLoaded = True
           
        
    def _get_ext_pep_chop_scores(self, ext_pep):
        
        tempDir = self._netChopTemp
        tempChopTxtIn  = os.path.join(tempDir, "%evd_temp_in.txt")
        tempChopTxtOut = os.path.join(tempDir, "%evd_temp_out.txt")
        tempProt = "Test prot"
        
        if not os.path.isdir(tempDir):
            os.mkdir(tempDir) 
            
        with open(tempChopTxtIn, "wb") as f:
            f.write(">%s\n%s"%(tempProt, ext_pep))   
        
        if not os.path.isfile(tempChopTxtOut):
            args = [self._binName, '-tdir', tempDir, tempChopTxtIn]
#             args = [self._binName, tempChopTxtIn]
            my_env = os.environ.copy()
            my_env["NETCHOP"] = self._netchopEnv
#             my_env["TMPDIR"] = tempDir
#             my_env["UNIX"] = "Darwin"
#             my_env["AR"] = "x86_64"
#             my_env["NMHOME"] = "/Users/ewaldvandyk/bin/netchop-3.1"
#             my_env["PWD"] = os.getcwd()
#             print my_env["PWD"]
            with open(tempChopTxtOut, "wb") as f:
                proc = subprocess.Popen(args, env=my_env, stdout=f)
                proc.wait()

        dataDict = netchop_raw_file_2_Dict(tempChopTxtOut, [tempProt])
        #raise IOError
        shutil.rmtree(tempDir, ignore_errors=False)

        return dataDict[tempProt]

    def _get_prot_subset_keys(self, pep_protName):
        protKeys = self._proteomeDict.keys()
        if pep_protName == None:
            return protKeys
        
        filt_protKeys = []
        for curr_key in protKeys:
            if pep_protName.upper() in curr_key.upper():
                filt_protKeys.append(curr_key)
        return filt_protKeys
    
    def _get_reduced_proteome(self, protSubSet, pep_start_pos, pepLen):
        reduced_proteome = dict()
        for protKey in protSubSet:
            curr_protLen = len(self._proteomeDict[protKey])
            if pep_start_pos == None:
                seqStart    = 0
                seqEnd      = curr_protLen
            else:
                seqStart = max(0, pep_start_pos - Peptide.aaExtendLeft-1)
                seqEnd   = min(curr_protLen, pep_start_pos + pepLen + Peptide.aaExtendRight+1)
            if seqEnd - seqStart >= pepLen:             
                reduced_proteome[protKey] = self._proteomeDict[protKey][seqStart:seqEnd]
        return reduced_proteome
    
    def get_score(self, \
                  test_pep      = None, \
                  germ_pep      = None, \
                  ext_pep       = None, \
                  pep_protName  = None, \
                  pep_start_pos = None):
        
        pepLen = len(test_pep)
        if not ext_pep == None:
            if not test_pep in ext_pep:
                raise ValueError("test_pep is not a substring of ext_pep")
            pepStart = ext_pep.index(test_pep)
            pepEnd = pepStart+pepLen-1 
            chopScores = self._get_ext_pep_chop_scores(ext_pep)
            return (None, None, "%f"%chopScores[pepEnd])
        self._load_proteome()
        protSubSet = self._get_prot_subset_keys(pep_protName)
        if not protSubSet:
            warnings.warn("Protein name not found in reference fasta file.", UserWarning)
            return ("", "", "")
        reduced_prot_Dict = self._get_reduced_proteome(protSubSet, pep_start_pos, pepLen)
        if not reduced_prot_Dict:
            warnings.warn("Protein in reference fasta is smaller than the test peptide to which it belongs.", UserWarning)
            return ("", "", "")
        if germ_pep == None:
            hamming_Dict = get_hamming_dict(reduced_prot_Dict, test_pep)
        else:
            hamming_Dict = get_hamming_dict(reduced_prot_Dict, germ_pep)
        (pep_ext_list, pep_ext_germ_list, pep_germ_list) = get_ext_pep_list_from_hamming(reduced_prot_Dict, hamming_Dict, test_pep)
        

        num_ambiguous_peps = len(pep_ext_list)
        ambChopAntScores = list()
        for ext_pep in pep_ext_list:
            pepStart = ext_pep.index(test_pep)
            pepEnd = pepStart+pepLen-1 
            chopScores = self._get_ext_pep_chop_scores(ext_pep)
            ambChopAntScores.append(chopScores[pepEnd])
            
        ambChopSelfScores = list()
        for (i, ext_pep) in enumerate(pep_ext_germ_list):
            pepStart = ext_pep.index(pep_germ_list[i])
            pepEnd = pepStart+pepLen-1 
            chopScores = self._get_ext_pep_chop_scores(ext_pep)
            ambChopSelfScores.append(chopScores[pepEnd])
        
        if      self._tieBreak == "MAX":
            antChopScore = max(ambChopAntScores)
            I = ambChopAntScores.index(antChopScore)
        elif    self._tieBreak == "MIN":
            antChopScore = min(ambChopAntScores)
            I = ambChopAntScores.index(antChopScore)
        elif    self._tieBreak == "MEAN":
            antChopScore = sum(ambChopAntScores) / num_ambiguous_peps
            I = 0
        else:
            antChopScore = max(ambChopAntScores)
            I = ambChopAntScores.index(antChopScore)
        
        selfChopScore = ambChopSelfScores[I]
        selfChopPep = pep_germ_list[I]

        print "__Antigen_peptide__"
        print test_pep
        print "__SelfRef_peptide__"
        print selfChopPep
        print "__Atigen_peptide_extended__"
        print pep_ext_list[I]
        print antChopScore
        print "__SelfRef_peptide_extended__"
        print pep_ext_germ_list[I]
        print selfChopScore
        print

        return (selfChopPep, "%f"%(selfChopScore), "%f"%(antChopScore))
        
            
        
def netchop_raw_file_2_Dict(tempChopTxt, protIDs):
    resDict = dict()
    
#         fasta_protSeqIDs, fasta_geneIDs = self._get_prot_keys_from_fasta()
    
    isHeader = False
    isData = False
    with open(tempChopTxt, "rb") as f:
        protI = -1
        protKey = None
        valueList = []
        for lineI, line in enumerate(f):
            if lineI % 1000 == 0:
                print lineI
            if Proteome.rawFields[0] in line and Proteome.rawFields[1] in line:
                isHeader = True
                protI += 1
                if not protKey == None:
                    resDict[protKey] = numpy.array(valueList)
#                     protKey = (fasta_protSeqIDs[protI], fasta_geneIDs[protI])
                protKey = protIDs[protI]
                if resDict.has_key(protKey):
                    raise LookupError("Protein ambiguity: Multiple proteins have the same name")
                valueList = []
                continue
            if isHeader and line[0] == "-":
                isHeader = False
                isData = True
                continue
            if isData and line[0] == "-":
                isData = False
                continue
            if isData:
                chopScore = float(line.split()[3])
                valueList.append(chopScore)
        resDict[protKey] = numpy.array(valueList)
    return resDict

def get_ext_pep_list_from_hamming(prot_Dict, hamming_Dict, test_pep):
    pepLen = len(test_pep)
    dict_keys = hamming_Dict.keys()
    min_hamming = min(hamming_Dict[dict_keys[0]])
    ext_pep_list = list()
    germ_ext_pep_list = list()
    germ_pep_list = list()
    for key in dict_keys:
        protLen = len(prot_Dict[key])
        minHammingI = [i for i, x in enumerate(hamming_Dict[key]) if x == min_hamming]
        for pepI in minHammingI:
            ext_start   = max(0, pepI - Peptide.aaExtendLeft-1)
            ext_end     = min(protLen, pepI + pepLen + Peptide.aaExtendRight + 1)
            curr_ext_pep = prot_Dict[key][ext_start:pepI] + test_pep + prot_Dict[key][(pepI+pepLen):ext_end]
            curr_germ_ext_pep = prot_Dict[key][ext_start:ext_end]
            curr_germ_pep = prot_Dict[key][pepI:(pepI+pepLen)]
            ext_pep_list.append(curr_ext_pep)
            germ_ext_pep_list.append(curr_germ_ext_pep)
            germ_pep_list.append(curr_germ_pep)
    return (ext_pep_list, germ_ext_pep_list, germ_pep_list)
        

def get_hamming_dict(protDict, test_pep):
    hammingDict = dict()
    minHamming = len(test_pep)
    dict_keys = protDict.keys()
    protMinHammings = [minHamming for _ in dict_keys]
    for (i, protKey) in enumerate(dict_keys):
        hammingDict[protKey] = hammingWindow(test_pep, protDict[protKey])
        protMinHammings[i] = min(hammingDict[protKey])
   
    minHamming = min(protMinHammings)
        
    for (i, protKey) in enumerate(dict_keys):
        if protMinHammings[i] > minHamming:
            del hammingDict[protKey]
    return hammingDict
    


def hammingWindow(subStr, fullStr):
    subStrLen = len(subStr)
    fullStrLen = len(fullStr)
    numPep = fullStrLen-subStrLen+1
    hammings = [0 for _ in range(numPep)]
    for i in range(numPep):
        hammings[i] = hamming(subStr, fullStr[i:(i+subStrLen)])
    return hammings

def hamming(s1, s2):
    """Calculate the Hamming distance between two strings"""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))
        
        
        
        
        
        
        
        
        
        
        
