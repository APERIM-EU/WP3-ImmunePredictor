import inspect
import os
# import subprocess
# import shutil
import csv
# import io_evd.zipbin
import io_evd.tsv
import numpy
import protbin.cls
import netchop.cls
import netmhcpan.cls
# import matplotlib
# import matplotlib.pyplot as plt
from __builtin__ import False


class Proteome:
    dataSubDir = "data"
    pepLenRange = [8,9,10,11,12]
    def __init__(self,
                 netchopEnv=None, \
                 netchopBin=None, \
                 netchopVer=None, \
                 netMHCpanBin=None, \
                 netMHCpanVer=None, \
                 dataDir=None, \
                 protName=None, \
                 mhcName=None, \
                 pepLen=None, \
                 chopTh=None, \
                 netMHCTh=None):
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Proteome.dataSubDir))
            
        self._initLocals()
            
        self.set_netchopEnv(netchopEnv) 
        self.set_netchopBin(netchopBin)
        self.set_netchopVer(netchopVer)
        self.set_netMHCpanBin(netMHCpanBin)
        self.set_netMHCpanVer(netMHCpanVer)
        self.set_dataDir(dataDir)
        self.set_protName(protName)
        self.set_mhcName(mhcName)
        self.set_pepLen(pepLen)
        self.set_chopTh(chopTh)
        self.set_netMHCTh(netMHCTh)
        
    def _initLocals(self):
        self._netchopEnv = None
        self._netchopBin = None
        self._netchopVer = None
        self._netMHCpanBin = None
        self._netMHCpanVer = None
        self._dataDir = None
        self._protName = None
        self._mhcName = None
        self._pepLen = None
        self._chopTh = None
        self._netMHCTh = None
        
        self._netchopEnvFlag = False
        self._netchopBinFlag = False
        self._netchopVerFlag = False
        self._netMHCpanBinFlag = False
        self._netMHCpanVerFlag = False
        self._dataDirFlag = False
        self._protNameFlag = False
        self._mhcNameFlag = False
        self._pepLenFlag = False
        self._chopThFlag = False
        self._netMHCThFlag = False
        
        self._allFlagsCheck = False
        
    def _update_inputs(self):
        if not  self._check_all_inputs():
            return
        self._filterDir = os.path.join(self._dataDir, 'filter')
        self._chopDir = os.path.join(self._filterDir, self._netchopVer)
        self._netMHCDir = os.path.join(self._chopDir, self._netMHCpanVer)
        self._protDir = os.path.join(self._netMHCDir, self._protName)
        self._paramDir = os.path.join(self._protDir, "chopTh_%.3f_netTh_%.3f"%(self._chopTh, self._netMHCTh))
        self._pepLenDir = os.path.join(self._paramDir, "%02d"%(self._pepLen))
        
        self._dataFile = os.path.join(self._pepLenDir, "%s.tsv"%(self._mhcName))
    
    def _check_all_inputs(self):
        allChecked = \
        self._netchopEnvFlag and \
        self._netchopBinFlag and self._netchopVerFlag and \
        self._netMHCpanBinFlag and self._netMHCpanVerFlag and \
        self._dataDirFlag and self._protNameFlag and \
        self._mhcNameFlag and self._pepLenFlag and \
        self._chopThFlag and self._netMHCThFlag
        
        self._allFlagsCheck = allChecked
        
        return allChecked
    
    def _get_netMHC_nM_th(self, netMHCDict):
        nMScores = numpy.array([])
        for protKey in netMHCDict.keys():
            nMScores = numpy.concatenate((nMScores, netMHCDict[protKey]), axis=0)
        nMScores = numpy.sort(nMScores, axis=0)
        self._nMScoresRanked = nMScores
        nM_cutI = round(self._netMHCTh*len(nMScores))
        nM_th = nMScores[nM_cutI]
#         plt.plot(nMScores[:round(self._netMHCTh*len(nMScores))])
#         plt.grid()
#         plt.show()
        return nM_th
    
    def _get_filtered_dict(self, protDict, chopDict, netMHCDict, nMTh):
        filtDict = dict()
        keyList = ["Peptide", "Mhc", "Proteome", self._netchopVer, \
                   "%s(nM)"%self._netMHCpanVer,  \
                   "%s(rank %% of proteome)"%self._netMHCpanVer, "ProteinName"]
        
        for i in range(len(keyList)):
            filtDict[keyList[i]] = []
        
        totNumProteins = float(len(self._nMScoresRanked))
        for protKey in protDict.keys():
            seqEnd = self._pepLen
            for seqStart in range(len(netMHCDict[protKey])):
                chopScore = chopDict[protKey][seqEnd-1]
                netMHCScore = netMHCDict[protKey][seqStart]
                if (chopScore >= self._chopTh) and \
                   (netMHCScore <= nMTh):
                    filtDict[keyList[0]].append(protDict[protKey][seqStart:seqEnd])
                    filtDict[keyList[1]].append(self._mhcName)
                    filtDict[keyList[2]].append(self._protName)
                    filtDict[keyList[3]].append(chopScore)
                    filtDict[keyList[4]].append(netMHCScore)
                    filtDict[keyList[5]].append(100*(numpy.where(self._nMScoresRanked <= netMHCScore)[0][-1]+1)/totNumProteins)
                    filtDict[keyList[6]].append(protKey) 
                seqEnd += 1
        return filtDict, keyList
            
    def set_netchopEnv(self, netchopEnv):
        if netchopEnv==None:
            self._netchopEnvFlag = False
            self._update_inputs()
            return
        if not os.path.isdir(netchopEnv):
            raise IOError("Netchop environment directory does not exist")
        self._netchopEnv = netchopEnv
        self._netchopEnvFlag = True
        self._update_inputs()
    
    def set_netchopBin(self, netchopBin):
        if netchopBin==None:
            self._netchopBinFlag = False
            self._update_inputs()
            return
        if not os.path.isfile(netchopBin):
            raise IOError("Netchop binary does not exist")
        self._netchopBin = netchopBin
        self._netchopBinFlag = True
        self._update_inputs()
        
    def set_netchopVer(self, netchopVer):
        if netchopVer==None:
            self._netchopEnvFlag = False
            self._update_inputs()
            return
        self._netchopVer = netchopVer
        self._netchopVerFlag = True
        self._update_inputs()
        
    def set_netMHCpanBin(self, netMHCpanBin):
        if netMHCpanBin==None:
            self._netMHCpanBinFlag = False
            self._update_inputs()
            return
        if not os.path.isfile(netMHCpanBin):
            raise IOError("NetMHCpan binary does not exist")
        self._netMHCpanBin = netMHCpanBin
        self._netMHCpanBinFlag = True
        self._update_inputs()
        
    def set_netMHCpanVer(self, netMHCpanVer):
        if netMHCpanVer==None:
            self._netMHCpanVerFlag = False
            self._update_inputs()
            return
        self._netMHCpanVer = netMHCpanVer
        self._netMHCpanVerFlag = True
        self._update_inputs()
        
    def set_dataDir(self, dataDir):   
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Proteome.dataSubDir))
            
        if not os.path.isdir(dataDir):
            if self._dataDir==None:
                self._dataDir = dataDir
                self._dataDirFlag = False
                self._update_inputs()
                return
            else:
                raise IOError("Data directory does not exist.")
        self._dataDir = dataDir
        self._dataDirFlag = True
        self._update_inputs()
            
    def set_protName(self, protName):
        if protName==None:
            self._protNameFlag = False
            self._update_inputs()
            return
        self._protName = protName
        self._protNameFlag = True
        self._update_inputs()
        
    def set_mhcName(self, mhcName):
        if mhcName==None:
            self._mhcNameFlag = False
            self._update_inputs()
            return
        if not (self._netMHCpanBinFlag and self._netMHCpanVerFlag and self._dataDirFlag):
            raise IOError("netMHCpan parameters must be set before selecting MHC name")
        mhcNnStream = netmhcpan.cls.Nn(binName=self._netMHCpanBin, \
                                             version=self._netMHCpanVer, \
                                             dataDir=self._dataDir)
        self._mhcName = mhcNnStream.get_NN(mhcName)
        self._mhcNameFlag = True
        self._update_inputs()
 
    def set_pepLen(self, pepLen):
        if pepLen==None:
            self._pepLenFlag = False
            self._update_inputs()
            return
        if not pepLen in Proteome.pepLenRange:
            raise IOError("Specified peptide length not in allowed range: %s"%(str(Proteome.pepLenRange)))
        self._pepLen = pepLen
        self._pepLenFlag = True
        self._update_inputs()
        
    def set_chopTh(self, chopTh):
        if chopTh==None:
            chopTh = 0.5
        else:
            chopTh = float(chopTh)
        if (chopTh < 0.0) or (chopTh > 1.0):
            raise IOError("chopTh has to be between 0.0 and 1.0")
        self._chopTh = chopTh
        self._chopThFlag = True
        self._update_inputs() 
        
    def set_netMHCTh(self, netMHCTh):
        if netMHCTh==None:
            netMHCTh = 0.1
        else:
            netMHCTh = float(netMHCTh)
        if (netMHCTh < 0.0) or (netMHCTh > 1.0):
            raise IOError("netMHCTh has to be between 0.0 and 1.0")
        self._netMHCTh = netMHCTh
        self._netMHCThFlag = True
        self._update_inputs()  
        
    def get_filtered_peps(self):
        if not self._allFlagsCheck:
            raise IOError("Not all input parameters set")
        if os.path.isfile(self._dataFile):
            tsvStream = io_evd.tsv.Data(inFile=self._dataFile)
            filtTsvDict = tsvStream.get_data()
            filtDict = io_evd.tsv.tsvDict2fieldDict(filtTsvDict)
            filtDict[self._netchopVer] = [float(i) for i in filtDict[self._netchopVer]]
            filtDict["%s(nM)"%self._netMHCpanVer] = [float(i) for i in filtDict["%s(nM)"%self._netMHCpanVer]]
            filtDict["%s(rank %% of proteome)"%self._netMHCpanVer] = [float(i) for i in filtDict["%s(rank %% of proteome)"%self._netMHCpanVer]]
            return filtDict
        if not os.path.isdir(self._pepLenDir):
            os.makedirs(self._pepLenDir)
        
        protBinStream = protbin.cls.Proteome(dataDir=self._dataDir, protName=self._protName)
        protDict = protBinStream.get_prots()
        
        chopBinStream = netchop.cls.Proteome(netchopEnv=self._netchopEnv, \
                                             binName=self._netchopBin, \
                                             version=self._netchopVer, \
                                             dataDir=self._dataDir, \
                                             protName=self._protName)
        chopDict = chopBinStream.get_scores()
        
        netMHCBinStream = netmhcpan.cls.Proteome(binName=self._netMHCpanBin, \
                                                 version=self._netMHCpanVer, \
                                                 dataDir=self._dataDir, \
                                                 protName=self._protName, \
                                                 mhcName=self._mhcName, \
                                                 pepLen=self._pepLen)
        netMHCDict = netMHCBinStream.get_scores()
        
        nMTh = self._get_netMHC_nM_th(netMHCDict)
        filtDict, filtKeyList = self._get_filtered_dict(protDict, chopDict, netMHCDict, nMTh)
        
        filtTsvDict = io_evd.tsv.fieldDict2tsvDict(filtDict)
#         print filtTsvDict
        
        tsvStream = io_evd.tsv.Data(outFile=self._dataFile)
        tsvStream.set_dialect(csv.excel_tab)
        tsvStream.set_data(filtTsvDict, keyList=filtKeyList)
#         tsvStream.set_data(filtTsvDict)
        tsvStream.write()
#         
        return filtDict

#         print filtDict.keys()
#         pepI = 100
#         print filtDict["ProteinName"][pepI]
#         print filtDict["Peptide"][pepI]
#         print filtDict["Cterm3_0"][pepI]
#         print filtDict["netMHCpan-3.0"][pepI]
#         plt.plot(filtDict["Cterm3_0"])
#         plt.grid()
#         plt.show()

#         print protDict
#         print chopDict
#         print netMHCDict
        
        