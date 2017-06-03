# BigRedButton
A pipeline to identify clusters of homologous RNA structure motifs from sequence data. The pipeline is optimised for execution on massively parallel computing infrastructure, although a local computing version will soon be developped. 
In summary, BigRedButton performs the following operations from a set of quesry sequences in .fasta format: 
1. Calculating base-pairing probabilities
2. All versus all pairwise comparisons with DotAligner 
3. Normalisation of resulting similarity matrix
4. OPTICS clustering and cluster extraction 
5. Multiple structural alignment of clustered sequences with mLocaRNA 

## DotAligner
At the core of the pipeline lies DotAligner, a pairwise sequence alignment heuristic which considers the ensemble of sub-optimal RNA base pair probabilities for both sequences. 
DotAligner is very efficient at segregating structurally distinct RNA families. 

## Benchmark data generation
java/GenerateRFAMsubsets.jar is a tool to sample RFAM Seed alignment stochastically with sequence constraints. Specifically, this tool will randomly yet exhaustively mine RFAM data to extract a minimal/maximal amount of sequences within user provided sequence identity ranges and sizes. 
