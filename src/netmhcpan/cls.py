import inspect
import os
import subprocess
# import shutil
import io_evd.zipbin
import io_evd.tsv
import csv
import numpy
import warnings
import datetime
# from _curses import version

class Peptide:
    def __init__(self, binName="/home/ewald/bin/netMHCpan-2.8", version="netMHCpan-2.8"):
        
        self._binName = binName
        self._version = version
        classFileName = inspect.getfile(self.__class__)
        classDir = os.path.dirname(classFileName)
        self._tempDir = classDir
        curr_time = datetime.datetime.now();
        curr_time_str = curr_time.strftime("%Y_%m_%d_%H_%M_%S_%f")
        self._pepInFile  = os.path.join(self._tempDir, 'pepInTemp_%s.txt'%curr_time_str)
        self._pepOutFile = os.path.join(self._tempDir, 'pepOutTemp_%s.txt'%curr_time_str)
         
        
    def get_score(self, peptide, mhcName):
        peptide = peptide.strip()
        mhcName = mhcName.strip()
        pepLen = len(peptide)
        
        with open(self._pepInFile, 'wb') as f:
            f.write("%s\n"%(peptide))
        
        
        args = [self._binName, "-v", "0", "-a", mhcName, "-s", \
                    "0", "-rth", "0.50", "-rlt", "2.00", \
                    "-l", "%d"%(pepLen), "-xls", "1", "-p", '1', \
                    "-f", self._pepInFile, "-xlsfile", self._pepOutFile]
        proc = subprocess.Popen(args)
        proc.wait()
        
        tsvStream = io_evd.tsv.Data(inFile=self._pepOutFile, firstLine=1)
        valueDict = tsvStream.get_data()
        
        aff_nM   = valueDict[0]["nM"]
        aff_Rank = valueDict[0]["Rank"]
        
        os.remove(self._pepInFile)
        os.remove(self._pepOutFile)
        
        nMField = "%s (nM)"%(self._version)
        rankField = "%s (rank %%)"%(self._version)
        
        returnDict = dict()
        returnDict[nMField] = aff_nM
        returnDict[rankField] = aff_Rank
        return returnDict

class Proteome:


    rawFields = ("Pos", "Peptide", "ID", "core", "1-log50k", "nM", "Rank", "Ave", "NB")
    dataSubDir = "data"
    protFastaSubDir = "protFasta"
    pepLenRange = [8,9,10,11,12]
    
    def __init__(self, binName="/home/ewald/bin/netMHCpan-2.8", version="netMHCpan-2.8",dataDir=None, protName=None, mhcName="HLA-A02:01", pepLen=9):
        
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Proteome.dataSubDir))
        self._binName = binName
        self._version  = version
        self._dataDir = dataDir
        self._protName = protName
        self._mhcName = mhcName
        self._pepLen = pepLen
        self._setup_paths()
    
    def _setup_paths(self, protName=None):
        if self._protName == None:
            raise TypeError("'protName' parameter not specified")
        self._check_mhcName_valid()
        self._check_pepLen_valid()
            
        self._protFastaDir = os.path.join(self._dataDir, Proteome.protFastaSubDir)
        self._protFastaFile = os.path.join(self._protFastaDir, "%s.fasta"%(self._protName))
        self._netMhcDir = os.path.join(self._dataDir, self._version)
        self._protDir = os.path.join(self._netMhcDir, self._protName)
        self._pepLenDir = os.path.join(self._protDir, "%02d"%(self._pepLen))
        self._dataFile = os.path.join(self._pepLenDir, "%s.dat"%(self._mhcName))

        if not (os.path.isfile(self._dataFile) or os.path.isfile(self._protFastaFile)):
            raise IOError("Proteome '%s' not recognized. Please ensure that the following file exists:\n%s"%(self._protName, self._protFastaFile))
        
        if not os.path.isdir(self._netMhcDir):
            os.mkdir(self._netMhcDir)
        if not os.path.isdir(self._protDir):
            os.mkdir(self._protDir)
        if not os.path.isdir(self._pepLenDir):
            os.mkdir(self._pepLenDir)
    
    def _check_mhcName_valid(self):
        mhcName_valid =  True
        if not mhcName_valid:
            raise TypeError("'mhcName' parameter is not valid")
    
    def _check_pepLen_valid(self):
        if self._pepLen in Proteome.pepLenRange:
            pepLenValid = True
        else:
            pepLenValid = False
        if not pepLenValid:
            raise TypeError("'pepLen' parameter is not valid")
    
    def get_scores(self):
        if not os.path.isfile(self._dataFile):
            dataDict = self.create_data_dict()
        else:
            inStream = io_evd.zipbin.Stream(self._dataFile)
            dataDict = inStream.read()
        return dataDict
    
    def create_data_dict(self):
        if not os.path.isfile(self._protFastaFile):
            raise IOError("Proteome '%s' not recognized. Please ensure that the following file exists:\n%s"%(self._protName, self._protFastaFile))
        
#         tempDir = os.path.join(self._pepLenDir, "%s_temp"%self._protName)
        curr_time = datetime.datetime.now();
        curr_time_str = curr_time.strftime("%Y_%m_%d_%H_%M_%S_%f")
        tempNetMhcTxt = os.path.join(self._pepLenDir, "%s_temp_%s.txt"%(self._protName, curr_time_str))
          
        if not os.path.isfile(tempNetMhcTxt):
            args = [self._binName, "-v", "0", "-a", self._mhcName, "-s", \
                    "0", "-rth", "0.50", "-rlt", "2.00", \
                    "-l", "%d"%(self._pepLen), "-xls", "1",  \
                    "-f", self._protFastaFile, "-xlsfile", tempNetMhcTxt]

#             my_env = os.environ.copy()
#             my_env["NETCHOP"] = self._netchopEnv
#             FNULL = open(os.devnull, 'w')
#             my_env = os.environ.copy()
#             my_env["PATH"] = "/usr/local/bin:" + my_env["PATH"]
#             proc = subprocess.Popen(args, env=my_env)
            proc = subprocess.Popen(args)
            proc.wait()

        dataDict = self._netmhcpan_raw_file_2_Dict(tempNetMhcTxt)
          
        outStream = io_evd.zipbin.Stream(self._dataFile)
        outStream.write(dataDict)

        os.remove(tempNetMhcTxt)
  
        return dataDict
    
    def _netmhcpan_raw_file_2_Dict(self, tempNetMhcTxt):
        resDict = dict()
        
        protIDs = self._get_prot_keys_from_fasta()
        
        isData = False
        with open(tempNetMhcTxt, "rb") as f:
            protI = -1
            protKey = None
            protID = None
            valueList = []
            protIs_seen = set()
            for lineI, line in enumerate(f):
                if lineI % 1000 == 0:
                    print lineI

                if Proteome.rawFields[0] in line and Proteome.rawFields[1] in line:
                    isData = True
                    continue
                if isData:
                    currFields = line.split()
                    currPos = currFields[0]
                    currNetProt= currFields[2]
                    nMScore = float(currFields[5])
                    if currPos == "0":
                        if not protKey == None:
                            resDict[protKey] = numpy.array(valueList)
                        protI += 1
                        protID_Is = set(self._net_protID_2_fasta_id_I(currNetProt, protIDs))
                        protID_Is = protID_Is - protIs_seen
                        if len(protID_Is) == 0:
                            raise IOError("Protein ID repeated in netMHC text output")
                        protID_I  = min(protID_Is)
                        protKey = protIDs[protID_I]
                        protIs_seen.add(protID_I)
                        if resDict.has_key(protKey):
                            raise LookupError("Protein ambiguity: Multiple proteins have the same name")
                        valueList = []
                    valueList.append(nMScore)  
            resDict[protKey] = numpy.array(valueList)
            for i in range(protI+1, len(protIDs)):
                protKey = protIDs[i]
                if resDict.has_key(protKey):
                    raise LookupError("Protein ambiguity: Multiple proteins have the same name")
                resDict[protKey] = numpy.array([])
        return resDict
    
    def _net_protID_2_fasta_id_I(self, netProtID, fastaIDs):
        fastaIs = []
        for i, fastaID in enumerate(fastaIDs):
            netProtIDparts = netProtID.split("_")
            currProtMatch = True
            for netProtIDpart in netProtIDparts:
                currProtMatch = currProtMatch and (netProtIDpart in fastaID)
            if currProtMatch:
                fastaIs.append(i)
        if len(fastaIs) == 0:
            raise IOError("Fasta ID not found for protein name %s" % (netProtID))
        if len(fastaIs) > 1:
            warnings.warn("multiple fasta IDs found for %s." % (netProtID), UserWarning)
        return fastaIs
    
    def _get_prot_keys_from_fasta(self):
        protIDs = []
        
        with open(self._protFastaFile, 'rb') as f:
            for line in f:
                if line[0] == ">":
                    protIDs.append(line[1:].strip())
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
   
   
class Nn:
    dataSubDir = "data"
    suppAlleleFile = "allHlaAlleles.txt"
    suppAlleleField = "hlas"
    suppAlleleList = None
    protName = "HIV-1_11706"
    mhcNNmapFile = "mhcNNmap.tsv"
    pepLen = 9
    mapFields = ["mhcName", "NNmhcName"]
    protMaxAbsDiffTol = 0.1
    
    def __init__(self, binName="/home/ewald/bin/netMHCpan-2.8", version="netMHCpan-2.8",dataDir=None, protName=None, pepLen=None):
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Nn.dataSubDir))
        
        if Nn.suppAlleleList == None:
            alleleFile = os.path.join(dataDir, Nn.suppAlleleFile)
            ioStream = io_evd.tsv.Data(inFile=alleleFile)
            ioStream.set_dialect(csv.excel_tab)
            Nn.suppAlleleList = ioStream.get_field(Nn.suppAlleleField)
            

        if protName == None:
            protName = Nn.protName
        if pepLen == None:
            pepLen = Nn.pepLen
            
        self._binName = binName
        self._version  = version
        self._dataDir = dataDir
        self._protName = protName
        self._pepLen = pepLen
        
        self._netMhcDir = os.path.join(self._dataDir, self._version)
        if not os.path.isdir(self._netMhcDir):
            os.mkdir(self._netMhcDir)
        
        self._mhcNNmapFile = os.path.join(self._netMhcDir, Nn.mhcNNmapFile)
        
    def _get_existing_NN(self):
        if not os.path.isfile(self._mhcNNmapFile):
            return None
        ioStream = io_evd.tsv.Data(inFile=self._mhcNNmapFile)
        ioStream.set_in_fields(Nn.mapFields)
        ioStream.overwrite_defined_outputs()
        for inFields in ioStream:
            if inFields[0] == self._mhcName:
                return inFields[1]
        return None
    
    def _get_all_NNs(self):
        if not os.path.isfile(self._mhcNNmapFile):
            return None
        ioStream = io_evd.tsv.Data(inFile=self._mhcNNmapFile)
        ioStream.set_in_fields(Nn.mapFields[1])
        ioStream.overwrite_defined_outputs()
        all_nns = []
        for nn in ioStream:
            all_nns.append(nn[0])
        all_nns = set(all_nns)
        return all_nns
    
    def _scores_are_equal(self, protDict1, protDict2):
        compKeys = protDict1.keys();
        maxDiff = 0
        for currKey in compKeys:
            currScoreDiff = protDict1[currKey] - protDict2[currKey]
            currMax = numpy.amax(numpy.absolute(currScoreDiff))
            if currMax > maxDiff:
                maxDiff = currMax
        
        print "max diff", maxDiff    
        if maxDiff < Nn.protMaxAbsDiffTol:
            return True
        return False
    
    def _get_new_NN(self):
        allNnNames = self._get_all_NNs()
        compProt = Proteome(binName=self._binName, version=self._version, \
                            dataDir=self._dataDir, protName=self._protName, \
                            pepLen=self._pepLen, mhcName=self._mhcName)
        compProtDict = compProt.get_scores()
        
        for currNnName in allNnNames:
            currProt = Proteome(binName=self._binName, version=self._version, \
                            dataDir=self._dataDir, protName=self._protName, \
                            pepLen=self._pepLen, mhcName=currNnName)
            currProtDict = currProt.get_scores()
            if self._scores_are_equal(compProtDict, currProtDict):
                return currNnName
            del currProt
        return self._mhcName
    
    def _append_map(self):
        ioStream = io_evd.tsv.Data(inFile=self._mhcNNmapFile, outFile=self._mhcNNmapFile)
        ioStream.set_out_fields(Nn.mapFields)
        ioStream.append_out_fields([self._mhcName, self._mhcNnName])
        ioStream.write()
        
        
    def get_NN(self, mhcName="HLA-A02:01"):

        if not mhcName in Nn.suppAlleleList:
            raise IOError("MHC allele name not supported")
        self._mhcName = mhcName
        newNnMapDict = dict()
        newNnMapDict[Nn.mapFields[0]] = self._mhcName
        newNnMapDict[Nn.mapFields[1]] = self._mhcName
        if not os.path.isfile(self._mhcNNmapFile):
            ioStream = io_evd.tsv.Data(outFile=self._mhcNNmapFile)
            ioStream.set_data([newNnMapDict])
            ioStream.write()
            self._mhcNnName = self._mhcName
            return self._mhcNnName
        self._mhcNnName = self._get_existing_NN()
        if self._mhcNnName == None:
            self._mhcNnName = self._get_new_NN()
            self._append_map()
        return self._mhcNnName

    def get_all_NN(self):
        nnList = set()
        for hla in Nn.suppAlleleList:
            nnHla = self.get_NN(hla)
            nnList.add(nnHla)
            print nnHla
            
        nnList = list(nnList)
        nnList.sort()
        return nnList
    
    def get_list_NN(self, hlaList):
        nnList = set()
        for hla in hlaList:
            nnHla = self.get_NN(hla)
            nnList.add(nnHla)
            print nnHla
            
        nnList = list(nnList)
        nnList.sort()
        return nnList
        
    
    def get_extracted_NN(self):
        print self._mhcNNmapFile
        mhcNNmapStream = io_evd.tsv.Data(inFile=self._mhcNNmapFile)
        nn_list = list(set(mhcNNmapStream.get_field(Nn.mapFields[1])))
        nn_list.sort()
        return nn_list
        
        
def raw_2_std_prot_id(rawId):
        tags = rawId.split("_")
        if len(tags) < 3:
            raise ValueError("Unique protein ID cannot be reconstructed from netMHCpan output")
        newID = tags[0] + "|" + tags[1] + "|"
        return newID
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        