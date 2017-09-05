import inspect
import os
import numpy
import warnings
import io_evd.tsv
import filter.cls


class SelfSim:
    dataSubDir = "data"
    pepLenRange = [9]
    
#     filter class parameters
    filterChopTh = 0.5
    filterNetMHCTh = 0.05
    
    pep_contact_range = [2,8]
    pep_noChange_range = [4,5]
    pep_left_range = [2,4]
    pep_right_range = [5,8]
    max_change_left_range = 1
    max_change_right_range = 1
    pmbec_cov_thresh = 0.05 
    
    
    def __init__(self,
                 netchopEnv=None, \
                 netchopBin=None, \
                 netchopVer=None, \
                 netMHCpanBin=None, \
                 netMHCpanVer=None, \
                 dataDir=None, \
                 protName=None, \
                 chopTh=None, \
                 netMHCTh=None):
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', SelfSim.dataSubDir))
            
        self._initLocals()
            
        self.set_netchopEnv(netchopEnv) 
        self.set_netchopBin(netchopBin)
        self.set_netchopVer(netchopVer)
        self.set_netMHCpanBin(netMHCpanBin)
        self.set_netMHCpanVer(netMHCpanVer)
        self.set_dataDir(dataDir)
        self.set_protName(protName)
        self.set_chopTh(chopTh)
        self.set_netMHCTh(netMHCTh)
        
        
#         self._load_PMBEC_matrix()
        
    def _load_PMBEC_matrix(self):
        pmbec_mat_file = os.path.join(self._calisSelfDir, 'PMBEC.txt')
        tsvStream = io_evd.tsv.Data(inFile=pmbec_mat_file)
        rawTsvDict = tsvStream.get_data()
        pmbecDict = io_evd.tsv.tsvDict2fieldDict(rawTsvDict)
        
        self._aas = pmbecDict["aa"]
        num_aas = len(self._aas)
        self._pmbec_mat = numpy.empty([num_aas, num_aas])
        for i in range(20):
            for j in range(20):
                self._pmbec_mat[i,j] = pmbecDict[self._aas[j]][i]
                
        self._pmbec_loaded = True
#         print self._get_pmbec_cov("E", "L")
        
    def _get_pmbec_cov(self, aa1, aa2):
        i1 = self._aas.index(aa1.upper())
        i2 = self._aas.index(aa2.upper())
        return self._pmbec_mat[i1, i2]  
    
    def _min_pep_pmbec_cov(self, pep1, pep2):
        if not len(pep1) == len(pep2):
            raise IOError("Peptide lenghts not equal for PMBEC covariance computation")
        pepLen = len(pep1)
        minCov = 1
        for i in range(pepLen):
            aa1 = pep1[i]
            aa2 = pep2[i]
            currCov = self._get_pmbec_cov(aa1, aa2)
            minCov = min(minCov, currCov)
        
        return minCov
        
    def _initLocals(self):
        self._netchopEnv = None
        self._netchopBin = None
        self._netchopVer = None
        self._netMHCpanBin = None
        self._netMHCpanVer = None
        self._dataDir = None
        self._protName = None
        self._chopTh = None
        self._netMHCTh = None
        
        self._netchopEnvFlag = False
        self._netchopBinFlag = False
        self._netchopVerFlag = False
        self._netMHCpanBinFlag = False
        self._netMHCpanVerFlag = False
        self._dataDirFlag = False
        self._protNameFlag = False
        self._chopThFlag = False
        self._netMHCThFlag = False
        
        self._allFlagsCheck = False
        
        self._calisSelfDir = None
        self._pmbec_loaded = False
        self._aas = None
        self._pmbec_mat = None
        self._mhcName = None
        
        
    def _update_inputs(self):
        if not  self._check_all_inputs():
            return
        self._calisSelfDir = os.path.join(self._dataDir, 'calisSelf')
        if not self._pmbec_loaded:
            self._load_PMBEC_matrix()
        
    def _check_all_inputs(self):
        allChecked = \
        self._netchopEnvFlag and \
        self._netchopBinFlag and self._netchopVerFlag and \
        self._netMHCpanBinFlag and self._netMHCpanVerFlag and \
        self._dataDirFlag and self._protNameFlag and \
        self._chopThFlag and self._netMHCThFlag
        
        self._allFlagsCheck = allChecked
        
        return allChecked
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
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', SelfSim.dataSubDir))
            
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
        
#     def set_mhcName(self, mhcName):
#         if mhcName==None:
#             self._mhcNameFlag = False
#             self._update_inputs()
#             return
#         if not (self._netMHCpanBinFlag and self._netMHCpanVerFlag and self._dataDirFlag):
#             raise IOError("netMHCpan parameters must be set before selecting MHC name")
#         mhcNnStream = netmhcpan.cls.Nn(binName=self._netMHCpanBin, \
#                                              version=self._netMHCpanVer, \
#                                              dataDir=self._dataDir)
#         self._mhcName = mhcNnStream.get_NN(mhcName)
#         self._mhcNameFlag = True
#         self._update_inputs()
 
#     def set_pepLen(self, pepLen):
#         if pepLen==None:
#             self._pepLenFlag = False
#             self._update_inputs()
#             return
#         if not pepLen in SelfSim.pepLenRange:
#             raise IOError("Specified peptide length not in allowed range: %s"%(str(SelfSim.pepLenRange)))
#         self._pepLen = pepLen
#         self._pepLenFlag = True
#         self._update_inputs()
        
    def set_chopTh(self, chopTh):
        if chopTh==None:
            self._chopThFlag = False
            self._update_inputs()
            return
        chopTh = float(chopTh)
        if (chopTh < 0.0) or (chopTh > 1.0):
            raise IOError("chopTh has to be between 0.0 and 1.0")
        self._chopTh = chopTh
        self._chopThFlag = True
        self._update_inputs() 
        
    def set_netMHCTh(self, netMHCTh):
        if netMHCTh==None:
            self._netMHCThFlag = False
            self._update_inputs()
            return
        netMHCTh = float(netMHCTh)
        if (netMHCTh < 0.0) or (netMHCTh > 1.0):
            raise IOError("netMHCTh has to be between 0.0 and 1.0")
        self._netMHCTh = netMHCTh
        self._netMHCThFlag = True
        self._update_inputs()  

    def get_degenerate_overlap_count(self, peptide, mhcName):
        pepLen = len(peptide)
        if not pepLen in SelfSim.pepLenRange:
            warnings.warn("Calis self similarity only supports 9-mer peptides. None returned", UserWarning)
            return float('nan')
        filterStream = filter.cls.Proteome()
        filterStream.set_netchopEnv(self._netchopEnv)
        filterStream.set_netchopBin(self._netchopBin)
        filterStream.set_netchopVer(self._netchopVer)
        filterStream.set_netMHCpanBin(self._netMHCpanBin)
        filterStream.set_netMHCpanVer(self._netMHCpanVer)
        filterStream.set_dataDir(self._dataDir)
        filterStream.set_protName(self._protName)
        filterStream.set_mhcName(mhcName)
        filterStream.set_pepLen(pepLen)
        filterStream.set_chopTh(SelfSim.filterChopTh)
        filterStream.set_netMHCTh(SelfSim.filterNetMHCTh)
        
        if self._mhcName == None or not self._mhcName == mhcName:
            self._filtDict = filterStream.get_filtered_peps()
            self._mhcName = mhcName
#         print self._filtDict.keys()
        numPeps = len(self._filtDict["Peptide"])
        
#         Filter unchanging aa
        domPepSymb = peptide[SelfSim.pep_noChange_range[0]:SelfSim.pep_noChange_range[1]]
        domAa = [self._filtDict["Peptide"][i][SelfSim.pep_noChange_range[0]:SelfSim.pep_noChange_range[1]]for i in range(numPeps)]
        domAaI = [i for i, x in enumerate(domAa) if x == domPepSymb]
        filtDict = make_filt_dict_with_index(self._filtDict, domAaI)
        
#         Filter netMHC rank
        netMHCThPerc = self._netMHCTh*100
        domAaI = [i for i, x in enumerate(filtDict["netMHCpan-3.0(rank % of proteome)"]) if x < netMHCThPerc]
        filtDict = make_filt_dict_with_index(filtDict, domAaI)
        
#         Filter left region
        numPeps = len(filtDict["Peptide"])
        domPepSymb = peptide[SelfSim.pep_left_range[0]:SelfSim.pep_left_range[1]]
        domAa = [filtDict["Peptide"][i][SelfSim.pep_left_range[0]:SelfSim.pep_left_range[1]] for i in range(numPeps)]
        domAaI = [i for i, x in enumerate(domAa) if hamdist(x, domPepSymb) <= SelfSim.max_change_left_range]
        filtDict = make_filt_dict_with_index(filtDict, domAaI)
        
#         Filter right region
        numPeps = len(filtDict["Peptide"])
        domPepSymb = peptide[SelfSim.pep_right_range[0]:SelfSim.pep_right_range[1]]
        domAa = [filtDict["Peptide"][i][SelfSim.pep_right_range[0]:SelfSim.pep_right_range[1]] for i in range(numPeps)]
        domAaI = [i for i, x in enumerate(domAa) if hamdist(x, domPepSymb) <= SelfSim.max_change_right_range]
        filtDict = make_filt_dict_with_index(filtDict, domAaI)
        
#         Filter min PBMEC cov
        numPeps = len(filtDict["Peptide"])
        domPepSymb = peptide[SelfSim.pep_contact_range[0]:SelfSim.pep_contact_range[1]]
        domAa = [filtDict["Peptide"][i][SelfSim.pep_contact_range[0]:SelfSim.pep_contact_range[1]] for i in range(numPeps)]
        domAaI = [i for i, x in enumerate(domAa) if self._min_pep_pmbec_cov(x, domPepSymb) >= SelfSim.pmbec_cov_thresh]
        filtDict = make_filt_dict_with_index(filtDict, domAaI)

#         Filter netChop
        domAaI = [i for i, x in enumerate(filtDict["Cterm3_0"]) if x >= self._chopTh]
        filtDict = make_filt_dict_with_index(filtDict, domAaI)
        
        uniPeps = set(filtDict["Peptide"])
        
        print len(uniPeps)
        return len(uniPeps)
        
def make_filt_dict_with_index(filtDict, I):
    new_dict = dict()
    new_dict["Peptide"] = [filtDict["Peptide"][i] for i in I]
    new_dict["netMHCpan-3.0(rank % of proteome)"] = [filtDict["netMHCpan-3.0(rank % of proteome)"][i] for i in I]
    new_dict["netMHCpan-3.0(nM)"] = [filtDict["netMHCpan-3.0(nM)"][i] for i in I]
    new_dict["Mhc"] = [filtDict["Mhc"][i] for i in I]
    new_dict["Cterm3_0"] = [filtDict["Cterm3_0"][i] for i in I]
    new_dict["ProteinName"] = [filtDict["ProteinName"][i] for i in I]
    new_dict["Proteome"] = [filtDict["Proteome"][i] for i in I]
        
    return new_dict
        
def hamdist(str1, str2):
    
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs        