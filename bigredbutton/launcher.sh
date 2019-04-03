#!/bin/bash

#################################################################
#                             Find binaries
#################################################################
# TO DO: Make more user friendy (i.e. check for presence)

module load gi/ViennaRNA/2.1.3 marsmi/dotaligner/0.3 marsmi/R/3.2.3

DOTALNBIN=$( which DotAligner ) 
RNAFOLDBIN=$( which RNAfold ) 
MLOCARNABIN=$( which mlocarna )
RSCRIPTBIN=$( which Rscript )

#################################################################
#################################################################
######
######             CUSTOM PARAMETERS BELOW --- EDIT ME
######
#################################################################
#################################################################

MAXGRIDCPU=179 # make an option for non HPC environments
GRID_ACCOUNT="RNABiologyandPlasticity"

#################################################################
#                            FUNCTIONS
#################################################################
>&2 echo -e "\e[91m             \`\`.:/osyyyys+:\`            "
>&2 echo -e "         \`.--:::::--::+oymNNms:         "
>&2 echo -e "       .-::-.-/sssyyyyyys++odNMm+\`      "
>&2 echo -e "     \`-::--odhsoosssyssooydmhhmNMm:     "
>&2 echo -e "    -:/-.smhyyhhhhhhhhhhhhhydNNNNMM+    "
>&2 echo -e "   :o+--mmhhhhddddddddddddhhhhmNNNMM+   "
>&2 echo -e "  .yy/.NmhdddddddddddddddddddddmNMNMN.  "
>&2 echo -e "  +dh-hNdddddddddd\e[39mBIG\e[91mdddddddddddmNMMMo  "
>&2 echo -e "  hNd/Nmdddddddddd\e[39mRED\e[91mdddddddddddhmMNMh  "
>&2 echo -e "  hNmyNmdddddddddd\e[39mBUTTON\e[91mdddddddhshMNNh  "
>&2 echo -e "  oMNmNNhdddddddddddddddddddddhh:mNNNo  "
>&2 echo -e "  .NMNNNdhhdddddddddddddddddhhho+MNNN.  "
>&2 echo -e "   +MMNMNdhhhdddddddddddddhhhh+/NNNN+   "
>&2 echo -e "    oMMNMNdhhhhhhdddddhhhhhho:yNNNN+    "
>&2 echo -e "     :mMNMMms+syhhhhhhhhyo//yNNNNm:     "
>&2 echo -e "      \`+mMMMMNds++////++ohNNNNNm+\`      "
>&2 echo -e "         :smMNNNNMMNNNNNNNNNds:         "
>&2 echo -e "            \`:+shddmmddhs+:\`            \e[39m"
>&2 echo -e                            


TIMEFORMAT="%R"
BASEDIR="$(cd "$(dirname "${1}")" && pwd)/$(basename "${1}")"
export BASEDIR=${BASEDIR%.*}

#################################################################
#                 Fold RNA to get dotplots
#################################################################
## TO DO : parallelize to make it faster 
if [[ ! -d $BASEDIR/logs ]]; then mkdir -p $BASEDIR/logs ; fi
#split fasta into smaller files
if [[ ! -d ${BASEDIR}/split ]]; then mkdir -p $BASEDIR/split ; fi
cd $BASEDIR/split
NUMSEQ=$(wc -l ${BASEDIR}/../$1 | awk '{print $1 }')
LINESperFILE=$(( ${NUMSEQ} / ${MAXGRIDCPU} ))
(( LINESperFILE+=1 ))

>&2 echo -ne "[*] Splitting fasta file with "$(( ${NUMSEQ} / 2 )) 
>&2 echo -ne " sequences into "${MAXGRIDCPU}" subfiles... "
split -l $LINESperFILE ${BASEDIR%/*}/$1
  >&2 echo -e "[ \e[92mDONE\e[39m ]"  

#prepare folding command
  CMD="qsub -V -terse -sync y -cwd -P ${GRID_ACCOUNT} -N fa2pp 
 -o $BASEDIR/logs/foldit."$(date +"%d%m%y-%H%M")".out 
 -e $BASEDIR/logs/foldit."$(date +"%d%m%y-%H%M")".err 
 -pe smp 1 
 -t 1:${MAXGRIDCPU}
 -b y -S /bin/bash
 ../fa2pp.sge"
cd $BASEDIR    && >&2 echo -e "[*] Working in BASEDIR: "$(pwd)

>&2 echo -ne "[*] Evaluating folder contents ..."
if [[ ! -d temp ]]; then mkdir temp ; fi
if [[ ! -d ps ]]; then mkdir ps ; fi
counter="-1"
if [[ ! -d pp ]]; then 
    mkdir pp
      >&2 echo -e "no RNA structures found"
      >&2 echo -e "    Generating pairwise probabilities..."
    $CMD | while read line ; do 
          counter=$(( $counter + 1 ))
          progress=$(( 100 * $counter / $MAXGRIDCPU ))
          echo -ne "\r... "$progress"% complete" 
      done
     >&2 echo -e "\r[*] Successfully completed command: "
     >&2 echo -e "\e[31m"$CMD"\e[39m" 
elif [[ $( ls pp | wc -l ) -ne $( grep ">" ../${1} | wc -l ) ]]; then 
    >&2 echo -e "incomplete amount of structures"
    cd $BASEDIR    
    >&2 echo -e "    Generating pairwise probabilities..."
    $CMD | while read line ; do 
          counter=$(( $counter + 1 ))
          progress=$(( 100 * $counter / $MAXGRIDCPU ))
          echo -ne "\r... "$progress"% complete" 
      done
     >&2 echo -e "\r[*] Successfully completed command: "
     >&2 echo -e "\e[31m"$CMD"\e[39m" 
else
    >&2 echo -e "found existing RNA structures "
      >&2 echo -ne "    \e[31mDo you want to refold them?\e[39m [y/N] "
    read -n1 redo 
    >&2 echo 
    if [[ $redo != 'n' && $redo != 'N' ]]; then 
         cd $BASEDIR    
         >&2 echo -e "[*] Generating pairwise probabilities..."
        $CMD | while read line ; do 
          counter=$(( $counter + 1 ))
          progress=$(( 100 * $counter / $MAXGRIDCPU ))
          echo -ne "\r... "$progress"% complete" 
      done
     >&2 echo -e "\r[*] Successfully completed command: "
     >&2 echo -e "\e[31m"$CMD"\e[39m" 
     fi
fi
cd $BASEDIR
rm -rf split
rm -rf temp/*

if [[ -e file_list.txt ]]; then rm file_list.txt ; fi
ls pp/ > file_list.txt && >&2 echo -e "[*] Saving list of files"


#################################################################
#                  Generate pairwise comparison list
#################################################################
>&2 echo -en "[*] Generating list of pairwise comparions (patience)...." 
cd ${BASEDIR}/pp
ls *pp | tee ../file_list.txt |\
  awk '
    OFS="\t" { arr[NR]=$0; ++numFiles } 
    END { for ( i=1; i <= numFiles; i++ ) 
      { for ( j=i; j <= numFiles ; j++) 
        print i,j,arr[i],arr[j]
      }
    }
  ' > $BASEDIR/pairwise_list.txt 
>&2 echo -e "... [ \x1B[92mDONE\x1B[39m ]"
cd $BASEDIR


#################################################################
#  SPLIT PAIRWISE LIST INTO $MAXCPU FILES FOR GRID COMPUTING
#################################################################
NUMPW=$(wc -l pairwise_list.txt | awk '{print $1}')
LINESperFILE=$(( ${NUMPW} / ${MAXGRIDCPU} ))
(( LINESperFILE+=1 ))
cd temp
# Clean up previous runs, if present
if [[ ! -d pairwise_groups ]]; then 
  mkdir pairwise_groups  
else 
  rm -rf pairwise_groups  && mkdir pairwise_groups
fi
cd pairwise_groups 
>&2 echo -n "[*] Splitting pairwise comparison list into "
>&2 echo -n $(( 1 + ${NUMPW}/$LINESperFILE ))" files of "
>&2 echo -n $LINESperFILE" lines" 
split -l $LINESperFILE $BASEDIR/pairwise_list.txt 
>&2 echo -e "... [ \x1B[92mDONE\x1B[39m ]"
cd $BASEDIR


################################################################
#                      RUN DOTALIGNER
################################################################
if [[ ! -d hpc ]]; then 
  mkdir hpc 
else 
  rm -rf hpc && mkdir hpc
fi
>&2 echo -e "[*] Performing all vs all pairwise comparisons "      
CMD="qsub -V -sync y -cwd -P ${GRID_ACCOUNT} -N DotAlnR -terse
 -o $BASEDIR/logs/launcher.$(date +"%d%m%y-%H%M").out 
 -e $BASEDIR/logs/launcher.$(date +"%d%m%y-%H%M").err 
 -pe smp 1 "
# Dotaligner doesnt need much memory
CMD=${CMD}" -l mem_requested=1256M,h_vmem=1512M
 -t 1:$((1+${NUMPW}/$LINESperFILE)) 
 -b y -S /bin/bash 
 ../worker.sge ${1}"
$CMD | while read line ; do 
      counter=$(( $counter + 1 ))
      progress=$(( 100 * $counter / $MAXGRIDCPU ))
      echo -ne "\r... "$progress"% complete" 
done

>&2 echo -e "\r[*] completed command: "
>&2 echo -e "\e[31m"$CMD"\e[39m" 



################################################################
#                 FILL IN DISTANCE MATRIX
################################################################
cd $BASEDIR
>&2 echo -n "[*] Collating scores..."
cd hpc 
cat *score > ../scores.txt && >&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

>&2 echo -n "[*] Concatenating alignment output..."
cat *out | gzip > ../alignments.out.gz 
>&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

>&2 echo -n "[*] Cleaning up intermediary files..."
cd .. 
rm -rf hpc && >&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

>&2 echo -n "[*] Normalising scores and converting to distance matrix..."
MAXMIN=$(  cat scores.txt | awk 'BEGIN{min = 2 ; max = -1} 
  { if ( $5 > max) {max = $5} else if ( $5 < min ) min = $5 } 
  END {print max" "min}'  )
MAX=${MAXMIN% *}
MIN=${MAXMIN#* }
awk -v min=$MIN -v max=$MAX 'OFS="\t"{ if ( $5 == "-0")  
  print $1,$2,$3,$4,"0" ; 
  else print $1,$2,$3,$4,($5-min)/(max-min) }' scores.txt  \
    > scores_normalized.tsv
awk 'OFS="\t"{ print $1,$2,$3,$4,1-$5 }'  scores_normalized.tsv \
    > dist.tsv
>&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

# get unique IDs; can make this prettier by stripping 
cut -f 3 dist.tsv | uniq > ids.txt 

################################################################
#                PREPARE CLUSTER ANALYSIS 
################################################################

>&2 echo -n "[*] Preparing cluster analysis script..."
cat > clustering.R << EOF
#!/usr/bin/Rscript
loadPackages <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , quiet = TRUE,
                        repos='http://cran.us.r-project.org' )
      #  Load package after installing
       suppressMessages( suppressWarnings(
                          require( i , character.only = TRUE )))
    }
  }
}
write("[R] loading required packages...", stderr())
suppressMessages( suppressWarnings(
    loadPackages( c("gplots","dbscan","RcolorBrewer","sparcl","data.table") ) ) ) 
write("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

write("[R] importing scores...", stderr())
# import dissimilarity scores
d <- read.table( "dist.tsv", header=F )
# this is faster
# d <- fread( "dist.tsv", header=F )
colnames( d ) <- c( 'x', 'y','x_name','y_name','score' )
write("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

write("[R] generating distance matrix...", stderr())
# convert to dissimilarity matrix
md <- matrix( nrow = max( d$y ), ncol=max( d$y ) )
md[ cbind( d$x , d$y ) ] <- as.numeric( d$score )
md[ cbind( d$y , d$x ) ] <- as.numeric( d$score )

# get unique names 
# d$Xname is of format "pp/PROTEIN_CELL_COORDINATES.pp"
# these next 5 lines pull out PROTEIN, which also 
# includes structured RNA controls and decoys
ids <- read.table("ids.txt", header=F)
#ids <- sapply(strsplit(as.character(ids),"/"), "[", 2)
#ids <- sapply(strsplit(as.character(ids),"\\."), "[", 1)

# convert to dist
rownames(md) <- ids$V1
colnames(md) <- ids$V1
D <- as.dist( md )
write("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

write("[R] performing density-based clustering ...", stderr())
#kNNdistplot(D,k=5)
Oc <- opticsXi( optics(D, eps=1, minPts=5, search="dist"), 
        xi = 0.006, minimum=T)
write("[ \x1B[92mDONE\x1B[39m ]\n", stderr())

write("[R] extracting clusters ...", stderr())
# print clusters to stdout
printClust <- function( O ) {
  for (cl in 1:NumClust) {
    print(paste("Cluster ======================== > ",cl)) ; 
    print(labels(D)[ O$cluster == cl ]) ;
  }
}
# write clusters out
dir.create( file.path( getwd(), "clusters" ), showWarnings = FALSE)
setwd("clusters")
NumClust <- max(Oc$cluster)
clusters <- 0 #make this 1 if optics(... , minimum=F )
for (cl in 1:NumClust) {
    l <- length(  Oc[Oc$cluster == cl] )
    # extract non-null clusters
    if ( l > 0  ) {
        sequences <- labels(D)[ Oc$cluster == cl ] 
        v <- sapply( strsplit( as.character( sequences ), "_"), "[", 1 ) 
        filename <- paste( "cluster_",cl,".txt", sep="" )
        write.table( sequences, file=filename , col.names=F, row.names=F, quote=F )
        clusters <- clusters + 1
    }
}
write("[ \x1B[92mDONE\x1B[39m ]\n", stderr())
write("Processed ",clusters," clusters\n")
write("[R] Processed **",clusters,"** clusters\n", stderr()) 
EOF
>&2 echo -e "[ \x1B[92mDONE\x1B[39m ]"

#################################################################
#                  Launch Cluster Analysis
#################################################################
>&2 echo -e "[*] Launching cluster analysis..."
Rscript --vanilla clustering.R
>&2 echo -e "[*] exiting :) " && exit 0 


#################################################################
#                  Process clusters from output
#################################################################
# extract fasta records from clusters
for file in cluster_*.txt ; do 
 mkdir ${file%*.txt}
 #change this to mv ? 
 cp $file ${file%*.txt}/
 cd ${file%*.txt}
 #TO DO filter out clusters with identical peaks 
 sed -e 's/_dp\.pp//g' -e 's/_/\.*/g' $file | while read line ; do echo '>'$line ; done > temp 
 grep -A 1 -f temp ${BASEDIR%/*}/$1 | grep -v -e '--' > ${file%*.txt}.fasta  
rm temp

$MLOCARNABIN --probabilistic --iterations=10 \
  --consistency-transformation --threads=6 --noLP ${file%*.txt}.fasta  
 cd ${file%*.txt}.out 
  ~/apps/ViennaRNA-2.2.10/bin/RNAalifold --color --aln -r results/result_prog.aln
  mv alirna.ps ${file%*.txt}_rna.ps
  mv aln.ps ${file%*.txt}_aln.ps
 cd ../..
 rm $file 
done
