# BigRedButton
All in one pipeline to identify clusters of homologous structures from a set on input sequences, such as CLIPseq peaks. Currently, the pipeline is designed to be launched on a high-performance computing infrastructure with Sun Grid Engine. 
A local version will soon be developped. 

## Usage

### From the command line: 
1. Copy bigredbutton folder contents to working directory
```
cd /path/to/workdir
cp /path/to/bigredbutton/bigredbutton/* ./
```

2. Make sure the launcher script is executable 
```
chmod 755 ./launcher.sh
```

3. Press the Big Red Button
```
./launcher.sh /path/to/sequences.fasta
```

### From the precompiled Docker image 
*In progress*
