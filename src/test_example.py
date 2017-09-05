import immGenPipe

# In toolParams you specify all the environmental variables required for netChop, netMHCpan and netMHCstabpan 
class toolParams:
#     Change these environmental variables for the specific system you use
#     NB! netChop does not work on OSX
#    netChopEnv = "/tbb/local/ii/src/netchop-3.1"
#    netChopBin = "/tbb/local/ii/src/netchop-3.1/bin/netChop"
#    netMHCpanBin = "/tbb/local/ii/src/netMHCpan-3.0/netMHCpan"
#    netMHCstabpanBin = "/tbb/local/ii/src/netMHCstabpan-1.0/netMHCstabpan"

    netChopEnv = "/Users/ewaldvandyk/bin/netchop-3.1/Darwin_x86_64"
    netChopBin = "/Users/ewaldvandyk/bin/netchop-3.1/Darwin_x86_64/bin/netChop"
    netMHCpanBin = "/tbb/local/ii/src/netMHCpan-3.0/netMHCpan"
    netMHCstabpanBin = "/tbb/local/ii/src/netMHCstabpan-1.0/netMHCstabpan"

#     Don't change these unless you really want to change the versions, In this case, 
#     all self proteome peptides will be re-extracted for hla alleles (time consuming)
    netChopVer = "Cterm3_0"
    netMHCpanVersion = "netMHCpan-3.0"
    netMHCstabpanVersion = "netMHCstabpan-1.0"

# Parameters used when calculating netChop scores for peptides
class netChopParams:
    #-------------------------------------------------------------------------
    ref_proteome = "HUMAN_9606"
                                    # Proteome fasta file name used to derive peptide context for netChop scores
                                    # For example, if the peptide is derived from HIV-1, an HIV-1 proteome should be provided
                                    # Proteome must be in FASTA format and located in <PIPE_DIR>/data/protFasta/
                                    # Note that the reference proteome should be representative of the antigen source, 
                                    # e.g. human for neoantigens and HIV for HIV derived peptides              
    #-------------------------------------------------------------------------
    tie_break  = "MAX"              
                                    # If it is unsure which protein is responsible for peptide, multiple netChop
                                    # scores will be computed. The reported netChop score is based on the rule
                                    # specified here. Possible options: {"MAX", "MIN", "MEAN"}
                                    # Default: "Max"

# Parameters used when computing self similarity scores for peptides.                              
class selfSimParams:
    #-------------------------------------------------------------------------
    self_proteome = "thymus_max60"  
                                    # Proteome fasta file name used for self similarity comparison. For example,
                                    # expressed proteins in the thymus
                                    # Proteome must be in FASTA format and located in <PIPE_DIR>/data/protFasta/
    chop_thresh   = None            
                                    # netChop threshold used when likely self peptides are selected. Default: 0.5
                                    # Set to "None" for default
    netMHC_thresh = None            
                                    # netMHCpan rank threshold used when likely self peptides are selected. Default: 0.023 (2.3%)
                                    # Set to "None" for default

# Input fields in the input file. Must be the same as the header names in the file. 
class inFields:
    #-------------------------------------------------------------------------
    mhcName         = "HLA"
                                    # Field with HLA allele name,e.g. "HLA-A02:01"
    #-------------------------------------------------------------------------
    test_pep        = "sequenceTum" 
                                    # Peptide for which scores are computed
    #-------------------------------------------------------------------------
    germ_pep        = None# "sequenceNorm"          
                                    # Field for reference germline peptide. Used when computing netChop scores for neoantigens.
                                    # For viruses where the germline peptide is the same, or if the field is missing, set germ_pep field to "None" 
                                    # "None" if field doesn't exist.                             
    #-------------------------------------------------------------------------
    ext_pep         = None
                                    # Field with extended aa sequence of test_pep. We recommend extending the peptide with at least 8 amino 
                                    # acids on both the N and C terminus side unless this extends beyond the protein start/end. The test peptide
                                    # must be a substring of the extended peptide. For neoantigens, this needs to come from the mutated version 
                                    # of the protein.
                                    # "None" if field doesn't exist
    #-------------------------------------------------------------------------
    pep_protName    = None #"pepProteinName"      
                                    # Field containing the name of the protein to which the test peptide belongs. 
                                    # "None" if field doesn't exist.
    #-------------------------------------------------------------------------
    pep_start_pos   = None          
                                    # Field indicating the start location of the germ_pep in the protein's reference FASTA sequence. 
                                    # Searches for closest matching peptide within +/-1 peptide start locations (to prevent off-by-one errors)
                                    # "None" to ignore the field.
    #-------------------------------------------------------------------------
                                #Note 1: "germ_pep", "ext_pep", "pep_protName" and "pep_start_pos" are only used when computing netChop scores
                                #        Although these fields are not strictly required, they help resolve protein ambiguity
                                #Note 2: if "ext_pep" is specified, "germ_pep", "pep_protName" and "pep_start_pos" will be ignored 
                                #Note 3: If germ_pep = None, test_pep will be used for germline matching
                                #Note 4: For multiple matches, ties are resolved according to the rule specified by "netChopParams.tie_break"
                                #Note 5: If germ_pep is specified and cannot be found in the reference proteome, netChop returns NA
                                #Note 6: If germ_pep is not provided, the germline peptide with the minimum hamming distance from "test_pep" 
                                #        is used for netChop scores. Ties are resolved according to the rule specified by "netChopParams.tie_break"
    


# Field to output
# netChop: netChop processing scores 
# netMHCpan: Compute netMHCpan affinity scores (and rank)
# netMHCstabpan: Compute netMHCstabpan stability scores (and rank)
# calisImmGen: Compute Calis Immunogenicity scores (PMID: 24204222). NB! Only 9-mer peptides supported, returns NaN otherwise
# calisSelf: Compute Calis self similarity scores (PMID: 22396638). NB! Only 9-mer peptides supported, returns NaN otherwise

class outFields:
    #Fields to output
    netChop         = True
    netMHCpan       = True
    netMHCstabpan   = True
    calisImmGen     = False
    calisSelf       = False
    
    
# Specify input and output file. Input file must be a tab/comma separated file with headers for each column
# Output file will have the same format with extra fields added 

#inFile = "/home/ewald/devel/python/immuGen/APERIM_WP3_pipeline_EVD/test/neoTon_test.txt"
#outFile = "/home/ewald/devel/python/immuGen/APERIM_WP3_pipeline_EVD/test/neoTon_out.txt"
inFile = "/Users/ewaldvandyk/devel/python/immGen/APERIM_WP3_pipeline_EVD/test/neoTon_test.txt"
outFile = "/Users/ewaldvandyk/devel/python/immGen/APERIM_WP3_pipeline_EVD/test/neoTon_out.txt"

featStream = immGenPipe.Feats(inFile=inFile, \
                              outFile=outFile, \
                              inFields=inFields, \
                              outFields=outFields, \
                              toolParams=toolParams, \
                              netChopParams=netChopParams, \
                              selfSimParams=selfSimParams)
featStream.add()
