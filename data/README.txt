All intermediate computations are stored in the data folder

All reference proteome fasta files should be placed in /data/protFasta/
All other subdirectories and data are handled automatically and should not be changed by the user

First time useage:
1. Download data.tar.gz from tbb.bio.uu.nl/ewald/immunePredict/
2. Unpack content into WP3_ImmunePredictor/data/ folder

All subsequent runs:
The data folder will naturally grow in size as more 'unseen' HLA alleles are processed in the pipeline
For a new 'unseen' HLA allele, predicting self similarity scores will be slow, since the peptidome needs to be predicted from the whole human thymus proteome. The newly extracted peptidome will stored in this data folder, which makes subesquent runs fast.
