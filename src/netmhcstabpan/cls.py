import inspect
import os
import subprocess
import io_evd.tsv
import datetime

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
                    "-1", "-rht", "0.5", "-rlt", "2.0", \
                    "-l", "%d"%(pepLen), "-xls", "1", "-p", '1', "-inptype", "1", \
                    "-f", self._pepInFile, "-xlsfile", self._pepOutFile]
        proc = subprocess.Popen(args)
        proc.wait()
        
        tsvStream = io_evd.tsv.Data(inFile=self._pepOutFile, firstLine=1)
        valueDict = tsvStream.get_data()

        stab_hl   = valueDict[0]["Thalf(h)"]
        stab_Rank = valueDict[0]["Rank"]
         
        os.remove(self._pepInFile)
        os.remove(self._pepOutFile)
         
        hlField = "%s (hours)"%(self._version)
        rankField = "%s (rank %%)"%(self._version)
         
        returnDict = dict()
        returnDict[hlField] = stab_hl
        returnDict[rankField] = stab_Rank
        return returnDict