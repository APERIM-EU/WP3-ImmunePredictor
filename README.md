ImmunePredictor is a pipeline developed in Python that predicts peptide cleavage, peptide-MHC affinity, stability, immunogenicity and self similarity (to human proteome) scores for human MHC class I.


The pipeline require the following external software tools:
1. NetChop 3.1 - Download: http://www.cbs.dtu.dk/services/NetChop/  - Note that we do not support the OSX version
2. NetMHCpan 3.0  - Download: http://www.cbs.dtu.dk/services/NetMHCpan/
3. NetMHCstabpan 1.0 - Download: http://www.cbs.dtu.dk/services/NetMHCstabpan/

Install and test:
1. Download the "WP3-ImmunePredictor/src/", "WP3-ImmunePredictor/data/" and "WP3-ImmunoPredictor/test/" directories. Make sure that "src", "data" and "test" are located in the same directory.
2. Download "data.tar.gz" from: http://tbb.bio.uu.nl/ewald/immunePredict/. Unpack the content into the "data" folder
3. Install NetChop, NetMHCpan and NetMHCstabpan
4. Open "WP3-ImmunePredictor/src/test_example.py"
5. In "test_example.py", locate the "toolParams" class and set the absolute paths for the netChop, netMHCpan and netMHCstabpan binaries. Note that netChop also requires an environment directory

Usage:
1. Run "test_example.py" to test the pipeline on an example input file from WP3-ImmunoPredictor/test/
2. See "test_example.py" for usage instructions
