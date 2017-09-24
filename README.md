# BigRedButton
A pipeline to identify clusters of homologous RNA structure motifs from sequence data. It can identify clusters of homologous RNA structure motifs from a set of user-provided sequences, such as RIPseq, CLIPseq, RNAse footprinting, etc. BigRedButton can identify known (and novel) RNA structures with good accuracy, as detailed in the paper. 

The pipeline is optimised for execution on massively parallel computing infrastructure, although a local computing version is being developped. In summary, BigRedButton performs the following operations from a set of query sequences in .fasta format: 
1. Calculating base-pairing probabilities
2. All versus all pairwise comparisons with DotAligner 
3. Normalisation of resulting similarity matrix
4. OPTICS clustering and cluster extraction 
5. Multiple structural alignment of clustered sequences with mLocaRNA 

The resulting multiple RNA structure alignments can then be used to calibrate a Covariance Model and scan a reference genome/transcriptome for homologs, thus improving the annotation of non-protein coding regions of the genome. 

More information on its usage, as well as installation of dependencies can be found in the ./bigredbutton directory. We also provide a Docker image to facilitate cross-platform operation. 
## DotAligner
At the core of the pipeline lies DotAligner, a fast pairwise sequence/structure alignment heuristic that considers the ensemble of sub-optimal RNA base pair probabilities for both sequences. DotAligner outperforms most more specialised RNA structure alignment algorithms, particularly within the 65-85% sequence identity range. DotAligner's can segregate structurally distinct RNA families with excellent accuracy.

More information on the usage of DotAligner for pairwise alignments can be found in the ./dotaligner directory.
## Benchmark data generation
To assess the performance of BigRedButton and other RNA structure clustering tools, we provide a tool to sample RFAM seed alignments stochastically, while allowing user-provided constraints on sequence composition, sequence length, and sample diversity.

The Java archive executable ./java/GenerateRFAMsubsets.jar will randomly yet exhaustively mine RFAM data to extract a minimal/maximal amount of sequences within user provided sequence identity ranges and sizes. 
