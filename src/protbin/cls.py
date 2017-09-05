import inspect
import os
import io_evd.zipbin


class Proteome:
    dataSubDir = "data"
    protFastaSubDir = "protFasta"
    protSubdirName = "protbin"
    def __init__(self, dataDir=None, protName=None):
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Proteome.dataSubDir))
        
        self._dataDir = dataDir
        self._setup_prot(protName)
    
    def _setup_prot(self, protName=None):
        if protName == None:
            raise TypeError("'protName' parameter not specified")
        self._protName = protName
        self._protFastaDir = os.path.join(self._dataDir, Proteome.protFastaSubDir)
        self._protFastaFile = os.path.join(self._protFastaDir, "%s.fasta"%(protName))
        self._protDataDir = os.path.join(self._dataDir, Proteome.protSubdirName)
        self._protDataFile = os.path.join(self._protDataDir, "%s.dat"%(self._protName))
        if not os.path.isdir(self._protDataDir):
            os.mkdir(self._protDataDir)
        
        if not (os.path.isfile(self._protDataFile) or os.path.isfile(self._protFastaFile)):
            raise IOError("Proteome '%s' not recognized. Please ensure that the following file exists:\n%s"%(self._protName, self._protFastaFile))
        
    def get_prots(self):
        if not os.path.isfile(self._protDataFile):
            dataDict = self.create_data_dict()
        else:
            inStream = io_evd.zipbin.Stream(self._protDataFile)
            dataDict = inStream.read()
        return dataDict
    
        
    def create_data_dict(self):
        if not os.path.isfile(self._protFastaFile):
            raise IOError("Proteome '%s' not recognized. Please ensure that the following file exists:\n%s"%(self._protName, self._protFastaFile))
         
        dataDict = self._fasta_2_dict()
        
        outStream = io_evd.zipbin.Stream(self._protDataFile)
        outStream.write(dataDict)
        
        return dataDict
    
    def _fasta_2_dict(self):
        resDict = dict()
        with open(self._protFastaFile, 'rb') as f:
            protI = -1
            protKey = None
            protString = ""
            for line in f:
                if line[0] == ">":
                    protI += 1
                    if not protKey == None:
                        resDict[protKey] = protString
#                     protSeqID = line[1:line.find("|",4)+1].strip()
#                     geneStart = line.find("GN=")+3
#                     geneEnd   = line.find("PE=")
#                     geneID    = line[geneStart:geneEnd].strip() 
#                     protKey = (protSeqID, geneID)
                    protKey = line[1:].strip()
                    if resDict.has_key(protKey):
                        raise LookupError("Protein ambiguity: Multiple proteins have the same name")
                    protString = ""
                else:
                    currString = line.strip(" \t\n").upper()
                    protString += currString
            resDict[protKey] = protString
        return resDict          
    
        
    def _fasta_2_dict_old(self):
        resDict = dict()
        with open(self._protFastaFile, 'rb') as f:
            protI = -1
            protKey = None
            protString = ""
            for line in f:
                if line[0] == ">":
                    protI += 1
                    if not protKey == None:
                        resDict[protKey] = protString
                    protSeqID = line[1:line.find("|",4)+1].strip()
                    geneStart = line.find("GN=")+3
                    geneEnd   = line.find("PE=")
                    geneID    = line[geneStart:geneEnd].strip()
                    protKey = (protSeqID, geneID)
                    if resDict.has_key(protKey):
                        raise LookupError("Protein ambiguity: Multiple proteins have the same name")
                    protString = ""
                else:
                    currString = line.strip(" \t\n").upper()
                    protString += currString
            resDict[protKey] = protString
        return resDict          
        

    