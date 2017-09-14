import inspect
import os
import subprocess
import io_evd.tsv
import datetime

class Peptide:
    dataSubDir = "data"
    
    def __init__(self, binName="/tbb/local/ii/src/netMHCstabpan-1.0/netMHCstabpan", version="netMHCstabpan-1.0", dataDir=None):
        if dataDir == None:
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            dataDir = os.path.abspath(os.path.join(classDir, '..', '..', Peptide.dataSubDir))
        self.set_bin_name(binName)
        self._version = version
        self._dataDir = dataDir
        self._stabpanDir = os.path.join(self._dataDir, self._version)
        if not os.path.isdir(self._stabpanDir):
            os.mkdir(self._stabpanDir)
        curr_time = datetime.datetime.now();
        curr_time_str = curr_time.strftime("%Y_%m_%d_%H_%M_%S_%f")
        self._pepInFile  = os.path.join(self._stabpanDir, 'pepInTemp_%s.txt'%curr_time_str)
        self._pepOutFile = os.path.join(self._stabpanDir, 'pepOutTemp_%s.txt'%curr_time_str)
         
    def set_bin_name(self, binName=None):
        if not os.path.isfile(binName):
            self._binName = None
            raise IOError("NetMHCstabpan bin file not accessible")
        else:
            self._binName = binName
        
    def get_score(self, peptide, mhcName):
        peptide = peptide.strip()
        mhcName = mhcName.strip()
        pepLen = len(peptide)
        
        with open(self._pepInFile, 'wb') as f:
            f.write("%s\n"%(peptide))
        
        try:
            args = [self._binName, "-v", "0", "-a", mhcName, "-s", \
                        "-1", "-rht", "0.5", "-rlt", "2.0", \
                        "-l", "%d"%(pepLen), "-xls", "1", "-p", '1', "-inptype", "1", \
                        "-f", self._pepInFile, "-xlsfile", self._pepOutFile]
            proc = subprocess.Popen(args)
            proc.wait()
            
            tsvStream = io_evd.tsv.Data(inFile=self._pepOutFile, firstLine=1)
            valueDict = tsvStream.get_data()
    
            stab_hl   = valueDict[0]["Thalf(h)"]
            stab_Rank = valueDict[0]["Rank"]
        except:
            os.remove(self._pepInFile)
            try:
                os.remove(self._pepOutFile)
            except:
                pass
            raise IOError("Unable to execute netMHCstabpan properly")
         
        os.remove(self._pepInFile)
        os.remove(self._pepOutFile)
         
        hlField = "%s (hours)"%(self._version)
        rankField = "%s (rank %%)"%(self._version)
         
        returnDict = dict()
        returnDict[hlField] = stab_hl
        returnDict[rankField] = stab_Rank
        return returnDict