# BigRedButton
All in one pipeline to identify clusters of homologous structures from a set on input sequences, such as CLIPseq peaks. Currently, the pipeline is designed to be launched on a high-performance computing infrastructure with Sun Grid Engine. 
A local multi-threaded version will soon be developped. 
N.B. for more than a few hundred sequences, a high-performance computing environment is recommended

## Usage

### From the command line: 
1. Copy bigredbutton folder contents to working directory
```
cd /path/to/workdir
cp /path/to/bigredbutton/bigredbutton/* ./
```
2. Edit the header lines of launcher.sh with the appropriate runtime parameters

3. Make sure the launcher script is executable 
```
chmod 755 ./launcher.sh
```

4. Press the Big Red Button
```
./launcher.sh /path/to/sequences.fasta
```

### From the precompiled Docker image 
*In progress*

## Under the hood
When using SGE, the pipeline will spawn array jobs with singe CPUs to (1) Fold input sequences and convert them to pairing probability matrices, and (2) perform pairwise sequence/structure alignments. These steps require little memory. 
For each cluster of RNA structure motifs detected, the pipeline then spawns a multiple structure alignment with mLorcaRNA using 6 CPUs, by default. 

## Output
Fasta files for each cluster will be written to the base directory where the pipeline was executed.
Multiple structure alignments (in Clustal .aln format) and visual representations of the structure motifs can be found within the sub-folders corresponding to each cluster. 
