import io_evd.tsv
import netchop.cls
import netmhcpan.cls
import netmhcstabpan.cls
import calis.immgen
import calis.cls
import time
# from __builtin__ import None

# toolParams:
# toolParams.netMHCpanBin -> full path + netMHCpan binary
# toolParams.netMHCpanVersion -> Version name, ex. "netMHCpan-2.8"
# toolParams.netMHCstabpanBin -> full path + netMHCstabpan binary
# toolParams.netMHCstabpanVersion -> Version name, ex. "netMHCstabpan-1.0"
#
# inFields:
# inFields.pepField -> Field name that contains peptide sequence
# inFields.mhcName -> Field name that contains MHC molecule names, ex "HLA-A02:01"


class inFields_default:
    def __init__(self):
        self.mhcName         = None
        self.test_pep        = None 
        self.ext_pep         = None
        self.germ_pep        = None          
        self.pep_protName    = None          
        self.pep_start_pos   = None          
 
class outFields_default:
    def __init__(self):
        self.netChop         = False
        self.netMHCpan       = False
        self.netMHCstabpan   = False
        self.calisImmGen     = False
        self.calisSelf       = False

class toolParams_default:
    def __init__(self):
        self.netChopEnv             = None
        self.netChopBin             = None
        self.netMHCpanBin           = None
        self.netMHCstabpanBin       = None
        self.netChopVer             = None
        self.netMHCpanVersion       = None
        self.netMHCstabpanVersion   = None
    
class netChopParams_default:
    def __init__(self):
        self.ref_proteome = None
        self.tie_break  = "MAX"
    
class selfSimParams_default:
    def __init__(self):
        self.self_proteome = None 
        self.chop_thresh   = 0.5            
        self.netMHC_thresh = 0.023            
 
class Feats:
    featOptions = ["netChop", "netMHCpan", "netMHCstabpan", "calisImmGen", "calisSelf"] 
    saveFrequency = 600#600 # Time in seconds
    
    def __init__(self, inFile=None, outFile=None, inFields=None, outFields=None, toolParams=None, netChopParams=None, selfSimParams=None):
        self._inFile = inFile
        self._outFile = outFile
        
        self.set_inFields(inFields)
        self.set_outFields(outFields)
        self.set_toolParams(toolParams)
        self.set_netChopParams(netChopParams)
        self.set_selfSimParams(selfSimParams)
                
    def _check_compulsory_params(self):  
        if self._inFields.mhcName == None:
            raise IOError("inFields.mhcName not specified")
        if self._inFields.test_pep == None:
            raise IOError("inFields.test_pep not specified")
        
        if self._toolParams.netChopEnv == None:
            raise IOError("toolParams.netChopEnv not specified")
        if self._toolParams.netChopBin == None:
            raise IOError("toolParams.netChopBin not specified")
        if self._toolParams.netMHCpanBin == None:
            raise IOError("toolParams.netMHCpanBin not specified")
        if self._toolParams.netMHCstabpanBin == None:
            raise IOError("toolParams.netMHCstabpanBin not specified")
        if self._toolParams.netChopVer == None:
            raise IOError("toolParams.netChopVer not specified")
        if self._toolParams.netMHCpanVersion == None:
            raise IOError("toolParams.netMHCpanVersion not specified")
        if self._toolParams.netMHCstabpanVersion == None:
            raise IOError("toolParams.netMHCstabpanVersion not specified")
        
        if (self._netChopParams.ref_proteome == None) and (self._outFields.netChop == True) and (self._inFields.ext_pep == None):
            raise IOError("netChopParams.ref_proteome not specified")
        
        if self._selfSimParams.self_proteome == None:
            raise IOError("selfSimParams.self_proteome not specified")
         
    def _get_inFields_list(self):
        self._inField_index = dict()
        self._inField_index["test_pep"] = 0;
        self._inField_index["mhc_Name"] = 1;
        index = 2
        
        inFields = [self._inFields.test_pep, self._inFields.mhcName]
        
        if self._inFields.germ_pep == None:
            self._inField_index["germ_pep"] = None
        else:
            inFields.append(self._inFields.germ_pep)
            self._inField_index["germ_pep"] = index
            index += 1    
        if self._inFields.ext_pep == None:
            self._inField_index["ext_pep"] = None
        else:
            inFields.append(self._inFields.ext_pep)
            self._inField_index["ext_pep"] = index
            index += 1
        if self._inFields.pep_protName == None:
            self._inField_index["pep_protName"] = None
        else:
            inFields.append(self._inFields.pep_protName)
            self._inField_index["pep_protName"] = index
            index += 1
        if self._inFields.pep_start_pos == None:
            self._inField_index["pep_start_pos"] = None
        else:
            inFields.append(self._inFields.pep_start_pos)
            self._inField_index["pep_start_pos"] = index
            index += 1
        return inFields   
               
    def _get_inFieldValues(self, element):
        inParams = inFields_default()
        if not self._inField_index["mhc_Name"] == None:
            inParams.mhcName = element[self._inField_index["mhc_Name"]]
        if not self._inField_index["test_pep"] == None:
            inParams.test_pep = element[self._inField_index["test_pep"]]
        if not self._inField_index["ext_pep"] == None:
            inParams.ext_pep = element[self._inField_index["ext_pep"]]
        if not self._inField_index["germ_pep"] == None:
            inParams.germ_pep = element[self._inField_index["germ_pep"]]
        if not self._inField_index["pep_protName"] == None:
            inParams.pep_protName = element[self._inField_index["pep_protName"]]
        if not self._inField_index["pep_start_pos"] == None:
            inParams.pep_start_pos = element[self._inField_index["pep_start_pos"]]
        return inParams    
             
    def set_inFields(self, inFields):
        self._inFields = inFields_default()
        if inFields == None:
            return
        
        if "mhcName" in dir(inFields) and not inFields.mhcName == None:
            self._inFields.mhcName = inFields.mhcName
        if "test_pep" in dir(inFields) and not inFields.test_pep == None:
            self._inFields.test_pep = inFields.test_pep
        if "ext_pep" in dir(inFields) and not inFields.ext_pep == None:
            self._inFields.ext_pep = inFields.ext_pep
        if "germ_pep" in dir(inFields) and not inFields.germ_pep == None:
            self._inFields.germ_pep = inFields.germ_pep
        if "pep_protName" in dir(inFields) and not inFields.pep_protName == None:
            self._inFields.pep_protName = inFields.pep_protName
        if "pep_start_pos" in dir(inFields) and not inFields.pep_start_pos == None:
            self._inFields.pep_start_pos = inFields.pep_start_pos
            
    def set_outFields(self, outFields):
        self._outFields = outFields_default()
        if outFields == None:
            return
        
        if "netChop" in dir(outFields) and not outFields.netChop == None:
            self._outFields.netChop = outFields.netChop
        if "netMHCpan" in dir(outFields) and not outFields.netMHCpan == None:
            self._outFields.netMHCpan = outFields.netMHCpan
        if "netMHCstabpan" in dir(outFields) and not outFields.netMHCstabpan == None:
            self._outFields.netMHCstabpan = outFields.netMHCstabpan
        if "calisImmGen" in dir(outFields) and not outFields.calisImmGen == None:
            self._outFields.calisImmGen = outFields.calisImmGen
        if "calisSelf" in dir(outFields) and not outFields.calisSelf == None:
            self._outFields.calisSelf = outFields.calisSelf
            
    def set_toolParams(self, toolParams):
        self._toolParams = toolParams_default()
        if toolParams == None:
            return
        
        if "netChopEnv" in dir(toolParams) and not toolParams.netChopEnv == None:
            self._toolParams.netChopEnv = toolParams.netChopEnv
        if "netChopBin" in dir(toolParams) and not toolParams.netChopBin == None:
            self._toolParams.netChopBin = toolParams.netChopBin
        if "netMHCpanBin" in dir(toolParams) and not toolParams.netMHCpanBin == None:
            self._toolParams.netMHCpanBin = toolParams.netMHCpanBin
        if "netMHCstabpanBin" in dir(toolParams) and not toolParams.netMHCstabpanBin == None:
            self._toolParams.netMHCstabpanBin = toolParams.netMHCstabpanBin
        if "netChopVer" in dir(toolParams) and not toolParams.netChopVer == None:
            self._toolParams.netChopVer = toolParams.netChopVer
        if "netMHCpanVersion" in dir(toolParams) and not toolParams.netMHCpanVersion == None:
            self._toolParams.netMHCpanVersion = toolParams.netMHCpanVersion
        if "netMHCstabpanVersion" in dir(toolParams) and not toolParams.netMHCstabpanVersion == None:
            self._toolParams.netMHCstabpanVersion = toolParams.netMHCstabpanVersion

    def set_netChopParams(self, netChopParams):
        self._netChopParams = netChopParams_default()
        if netChopParams == None:
            return
        
        if "ref_proteome" in dir(netChopParams) and not netChopParams.ref_proteome == None:
            self._netChopParams.ref_proteome = netChopParams.ref_proteome
        if "tie_break" in dir(netChopParams) and not netChopParams.tie_break == None:
            self._netChopParams.tie_break = netChopParams.tie_break
            
    def set_selfSimParams(self, selfSimParams):
        self._selfSimParams = selfSimParams_default()
        if selfSimParams == None:
            return
        
        if "self_proteome" in dir(selfSimParams) and not selfSimParams.self_proteome == None:
            self._selfSimParams.self_proteome = selfSimParams.self_proteome
        if "chop_thresh" in dir(selfSimParams) and not selfSimParams.chop_thresh == None:
            self._selfSimParams.chop_thresh = selfSimParams.chop_thresh
        if "netMHC_thresh" in dir(selfSimParams) and not selfSimParams.netMHC_thresh == None:
            self._selfSimParams.netMHC_thresh = selfSimParams.netMHC_thresh
            
       
    def add(self):
        self._check_compulsory_params()
        inFields = self._get_inFields_list()
        
        tsvIoStream = io_evd.tsv.Data(inFile=self._inFile,outFile=self._outFile)
        tsvIoStream.set_in_fields(inFields)
        
        self.save_time_passed()
        if self._outFields.netChop:
            netChopStream = netchop.cls.Peptide(netchopEnv      = self._toolParams.netChopEnv, \
                                                binName         = self._toolParams.netChopBin, \
                                                version         = self._toolParams.netChopVer, \
                                                proteomeName    = self._netChopParams.ref_proteome, \
                                                tieBreak        = self._netChopParams.tie_break)
            
            refSelfPep_field    = "%s_refSelf_pep"%(self._toolParams.netChopVer)
            refSelfScore_field  = "%s_refSelf_score"%(self._toolParams.netChopVer)
            AntigenScore_field  = "%s_antigen_score"%(self._toolParams.netChopVer)
            
            if not self._inFields.ext_pep == None:
                tsvIoStream.set_out_fields([AntigenScore_field])
            else:
                tsvIoStream.set_out_fields([refSelfPep_field, refSelfScore_field, AntigenScore_field])
            for element in tsvIoStream:
                inParams = self._get_inFieldValues(element)
                (selfPep, selfScore, antScore) = netChopStream.get_score(test_pep = inParams.test_pep, \
                                                    germ_pep        = inParams.germ_pep, \
                                                    ext_pep         = inParams.ext_pep, \
                                                    pep_protName    = inParams.pep_protName, \
                                                    pep_start_pos   = inParams.pep_start_pos)
                if not self._inFields.ext_pep == None:
                    tsvIoStream.write_out_fields([antScore])
                else:
                    tsvIoStream.write_out_fields([selfPep, selfScore, antScore])
                if self.save_time_passed():
                    tsvIoStream.write()
                    
                
                    
        
        if self._outFields.netMHCpan:
            netMHCpanStream = netmhcpan.cls.Peptide(binName=self._toolParams.netMHCpanBin, version=self._toolParams.netMHCpanVersion)
            nM_field = "%s (nM)"%(self._toolParams.netMHCpanVersion)
            rankField = "%s (rank %%)"%(self._toolParams.netMHCpanVersion)
            tsvIoStream.set_out_fields([nM_field, rankField])
            for element in tsvIoStream:
                inParams = self._get_inFieldValues(element)
                affDict = netMHCpanStream.get_score(peptide=inParams.test_pep, mhcName=inParams.mhcName)
                tsvIoStream.write_out_fields([affDict[nM_field], affDict[rankField]])
                if self.save_time_passed():
                    tsvIoStream.write()
                
        if self._outFields.netMHCstabpan:
            netMHCstabpanStream = netmhcstabpan.cls.Peptide(binName=self._toolParams.netMHCstabpanBin, version=self._toolParams.netMHCstabpanVersion)
            hl_field = "%s (hours)"%(self._toolParams.netMHCstabpanVersion)
            rankField = "%s (rank %%)"%(self._toolParams.netMHCstabpanVersion)
            tsvIoStream.set_out_fields([hl_field, rankField])
            for element in tsvIoStream:
                inParams = self._get_inFieldValues(element)
                stabDict = netMHCstabpanStream.get_score(peptide=inParams.test_pep, mhcName=inParams.mhcName)
                tsvIoStream.write_out_fields([stabDict[hl_field], stabDict[rankField]])
                if self.save_time_passed():
                    tsvIoStream.write()
    
        if self._outFields.calisImmGen:
            immGen_field = "CalisImmGen"
            tsvIoStream.set_out_fields(immGen_field)
            for element in tsvIoStream:
                inParams = self._get_inFieldValues(element)
                immGenScore = calis.immgen.score(pepSeq=inParams.test_pep, HLA=inParams.mhcName)
                tsvIoStream.write_out_fields(str(immGenScore))
                if self.save_time_passed():
                    tsvIoStream.write()
                
        if self._outFields.calisSelf:
            self_field = 'Calis_self_degenerate_count'
            calisSelfSimStream = calis.cls.SelfSim()
            calisSelfSimStream.set_netchopEnv(self._toolParams.netChopEnv)
            calisSelfSimStream.set_netchopBin(self._toolParams.netChopBin)
            calisSelfSimStream.set_netMHCpanBin(self._toolParams.netMHCpanBin)
            calisSelfSimStream.set_netchopVer(self._toolParams.netChopVer)
            calisSelfSimStream.set_netMHCpanVer(self._toolParams.netMHCpanVersion)
            calisSelfSimStream.set_protName(self._selfSimParams.self_proteome)
            calisSelfSimStream.set_chopTh(self._selfSimParams.chop_thresh)
            calisSelfSimStream.set_netMHCTh(self._selfSimParams.netMHC_thresh)
            tsvIoStream.set_out_fields(self_field)
            for element in tsvIoStream:
                inParams = self._get_inFieldValues(element)
                degenScore = calisSelfSimStream.get_degenerate_overlap_count(peptide=inParams.test_pep, mhcName=inParams.mhcName)
                tsvIoStream.write_out_fields(str(degenScore))
                if self.save_time_passed():
                    tsvIoStream.write()
                
        tsvIoStream.write()
        
    def save_time_passed(self):
        try:
            self._prevSaveTime
        except AttributeError:
            self._prevSaveTime = time.time()
            return False
        currTime = time.time()
        if currTime - self._prevSaveTime > Feats.saveFrequency:
            self._prevSaveTime = currTime
            return True
        else:
            return False 
        
