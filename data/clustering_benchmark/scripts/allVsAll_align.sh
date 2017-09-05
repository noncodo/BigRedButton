#!/bin/bash
MAXGRIDCPU=179
TIMEFORMAT="%R"
GRID_ACCOUNT="RNABiologyandPlasticity"

# load modules (for or local server)
module load marsmi/locarna/1.7.13 gi/ViennaRNA/2.1.3 marsmi/dotaligner/0.3 marsmi/foldalign/2.1.1

# >&2 echo "Generating list of pairwise comparions...." 
# k=0
# Files=$( ls ${1}/pp/*.pp | tee ${1}/file_list.txt | wc -l )
# for (( i=1 ; i <= ${Files} ; i++ )); do 
#   for (( j = $i ; j <= ${Files}; j++ )); do 
#       F1=$( sed $i'q;d' ${1}/file_list.txt )
#       F2=$( sed $j'q;d' ${1}/file_list.txt )
#       echo -e $i"\t"$j"\t"`pwd`/${F1//pp/ps}"\t"`pwd`/${F2//pp/ps} 
#       ((k+=1))   
#       if [[ $(( $k % 50 )) -eq 0 ]]; then 
#         >&2 printf "\r%.1f%s" $( echo $k $Files | \
#           awk '{printf 100*$1/($2*($2+1)/2)}' ) "% of alignments completed"
#       fi
#   done
# done > $1/pairwise_list.txt
# >&2 echo  
# cut -d "/" -f 3 $1/file_list.txt | cut -d "_" -f 1-2 > $1/$1_structure_ids.txt

## # #############################################
## # #  MAKE BINARY MATRIX FOR CORRELATION
## # #############################################
##
##    N.B. this has been placed in the downstream R code instead

## >&2  echo -n "[NOTE] generating matrix..."
## LINES=$( wc -l $1/DA.log )
## k=1 
## ###   Print out column names 
## (sed 's/^.*\///g' low_pi/file_list.txt | cut -d "_" -f 1,2 | while read line ; 
##   do  echo -n $line" " 
##   done ; echo ) | sed 's/.$//' > $1/binary.matrix
#
## for (( i=1 ; i <= $Files ; i++ )); do 
##   S1=$(sed $i'q;d' $1/file_list.txt | sed 's/^.*\///' | cut -d "_" -f 1 )     
##   for (( j=1 ; j <= $Files; j++ )); do 
##     S2=$(sed $j'q;d' $1/file_list.txt | sed 's/^.*\///' | cut -d "_" -f 1 ) 
##     ###   Print out row names 
##     #if [[ $j -eq 1 ]]; then 
##     #    echo -n $S1" " >> $1/binary.matrix
##     if [ $S1 == $S2 ]; then 
##       echo -n "1 " 
##     else
##       echo -n "0 " 
##     fi
##   done
##   echo 
## done | sed 's/.$//' >>  $1/binary.matrix
## >&2 echo "DONE"

#################################################################
# SPLIT PAIRWISE LIST INTO $MAXCPU FILES FOR GRID COMPUTING
#################################################################
NUMPW=$(wc -l $1/pairwise_list.txt | awk '{print $1}')
LINESperFILE=$(( ${NUMPW} / ${MAXGRIDCPU} ))
(( LINESperFILE+=1 ))
# cd $1
# if [[ ! -d pairwise_groups ]]; then 
#   mkdir pairwise_groups  
# else rm pairwise_groups/* 
# fi
#  cd pairwise_groups 
#  >&2 echo "Splitting pairwise comparison list into "$((1+${NUMPW}/$LINESperFILE))" files of "$LINESperFILE" lines" \
#   && split -l $LINESperFILE ../pairwise_list.txt 
#  cd ..
# cd ..

# if [[ ! -d $1/logs ]]; then mkdir $1/logs ; fi
#################################################################
#                       RUN DOTALIGNER
#################################################################
##   N. B. 
##   This was for initial benchmarking. 
##   A subsequent fine-tuned parameter optimisation was done using rfam_param_opt

 # if [[ ! -d $1/dotaligner ]]; then mkdir $1/dotaligner ; fi      # pre-existing files get removed in the sub-commands
 # if [[ -e ${1}/qsub.dotalnr.out ]] ; then rm ${1}/qsub.dotalnr.out ; fi 
 # if [[ -e ${1}/qsub.dotalnr.err ]] ; then rm ${1}/qsub.dotalnr.err ; fi 
 # CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N DotAlnR -o ${1}/logs/dotaligner.$(date +"%d%m%y-%H%M").out -e ${1}/logs/dotaligner.$(date +"%d%m%y-%H%M").err -pe smp 1"
 # CMD=${CMD}" -l mem_requested=1256M,h_vmem=1512M,h=!(epsilon-0-24.local|epsilon-0-25.local) -t 1:$((1+${NUMPW}/$LINESperFILE)) -b y ./dotaligner.sge ${1}"
 # DOTALNR=$( $CMD ) 
 # echo "launching: "$CMD 
 # echo $DOTALNR
 # JID=$( echo ${DOTALNR} | cut -d " " -f 3 | cut -d "." -f 1 )
 # qsub -hold_jid $JID -V -o ./$1/logs/dotaligner.$(date +"%d%m%y-%H%M").clean -cwd -N clean_dotalnr -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 dotaligner


#################################################################
#                         RUN LOCARNA
#################################################################
# if [[ ! -d $1/locarna ]]; then mkdir $1/locarna ; fi            # pre-existing files get removed in the sub-commands
# if [[ -e ${1}/qsub.locarna.out ]] ; then rm ${1}/qsub.locarna.out ; fi 
# if [[ -e ${1}/qsub.locarna.err ]] ; then rm ${1}/qsub.locarna.err ; fi 
# CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N LOCARNA -o ${1}/logs/locarna.$(date +"%d%m%y-%H%M").out -e ${1}/logs/locarna.$(date +"%d%m%y-%H%M").err -pe smp 1"
# CMD=${CMD}" -l mem_requested=1792M,h_vmem=2048M -t 1:$((1+${NUMPW}/$LINESperFILE)) -b y ./locarna.sge ${1}"
# LOCARNA=$( $CMD ) 
# echo "launching: "$CMD 
# echo $LOCARNA
# JID=$( echo ${LOCARNA} | cut -d " " -f 3 | cut -d "." -f 1 )
# qsub -hold_jid $JID -V -o ./$1/logs/locarna.$(date +"%d%m%y-%H%M").clean -cwd -N clean_locarna -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 locarna


#################################################################
#                         RUN CARNA
#################################################################
 # if [[ ! -d $1/carna ]]; then mkdir $1/carna ; fi                 # pre-existing files get removed in the sub-commands
 # if [[ -e ${1}/qsub.carna.out ]] ; then rm ${1}/qsub.carna.out ; fi 
 # if [[ -e ${1}/qsub.carna.err ]] ; then rm ${1}/qsub.carna.err ; fi 
 # CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N CARNA -o ${1}/logs/carna.$(date +"%d%m%y-%H%M").out -e ${1}/logs/carna.$(date +"%d%m%y-%H%M").err -pe smp 1"
 # CMD=${CMD}" -l  mem_requested=2792M,h_vmem=3048M -t 1:$((1+${NUMPW}/${LINESperFILE})) -b y ./carna.sge ${1}"
 # CARNA=$( $CMD ) 
 # echo "launching: "$CMD
 # echo $CARNA
 # #Example output: 
 # #Your job-array 4425908.1-192:1 ("CARNA") has been submitted
 # # N.B. qsub option "-terse" only outputs the job number
 # JID=$( echo ${carnaA} | cut -d " " -f 3 | cut -d "." -f 1 )
 # qsub -hold_jid $JID -V -o ./$1/logs/carna.$(date +"%d%m%y-%H%M").clean -cwd -N clean_carna -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 carna


#################################################################
#                         RUN PMCOMP
#################################################################
 # if [[ ! -d $1/pmcomp ]]; then mkdir $1/pmcomp ; fi                 # pre-existing files get removed in the sub-commands
 # if [[ -e ${1}/qsub.pmcomp.out ]] ; then rm ${1}/qsub.pmcomp.out ; fi 
 # if [[ -e ${1}/qsub.pmcomp.err ]] ; then rm ${1}/qsub.pmcomp.err ; fi 
 # CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N pmcomp -o ${1}/logs/pmcomp.$(date +"%d%m%y-%H%M").out -e ${1}/logs/pmcomp.$(date +"%d%m%y-%H%M").err -pe smp 1"
 # CMD=${CMD}" -l  mem_requested=7G,h_vmem=7G -t 1:$((1+${NUMPW}/${LINESperFILE})) -b y ./pmcomp.sge ${1}"
 # pmcomp=$( $CMD ) 
 # echo "launching: "$CMD
 # echo $pmcomp
 # #Example output: 
 # #Your job-array 4425908.1-192:1 ("pmcomp") has been submitted
 # # N.B. qsub option "-terse" only outputs the job number
 # JID=$( echo ${pmcomp} | cut -d " " -f 3 | cut -d "." -f 1 )
 # qsub -hold_jid $JID -V -o ./$1/logs/pmcomp.$(date +"%d%m%y-%H%M").clean -cwd -N clean_pmcomp -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup_pmcomp.sge $1 pmcomp


################################################################
#                         RUN Foldalign
#################################################################
 if [[ ! -d $1/foldalign ]]; then mkdir $1/foldalign ; fi                 # pre-existing files get removed in the sub-commands
 if [[ -e ${1}/qsub.foldalign.out ]] ; then rm ${1}/qsub.foldalign.out ; fi 
 if [[ -e ${1}/qsub.foldalign.err ]] ; then rm ${1}/qsub.foldalign.err ; fi 
 CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N foldalign -o ${1}/logs/foldalign.$(date +"%d%m%y-%H%M").out -e ${1}/logs/foldalign.$(date +"%d%m%y-%H%M").err -pe smp 1"
 CMD=${CMD}" -l  mem_requested=2.5G,h_vmem=2.5G -t 1:$((1+${NUMPW}/${LINESperFILE})) -b y ./foldalign.sge ${1}"
 foldalign=$( $CMD ) 
 echo "launching: "$CMD
 echo $foldalign
 #Example output: 
 #Your job-array 4425908.1-192:1 ("foldalign") has been submitted
 # N.B. qsub option "-terse" only outputs the job number
 JID=$( echo ${foldalign} | cut -d " " -f 3 | cut -d "." -f 1 )
 qsub -hold_jid $JID -V -o ./$1/logs/foldalign.$(date +"%d%m%y-%H%M").clean -cwd -N clean_foldalign -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 foldalign

################################################################
#                         RUN RNApaln
################################################################
 # if [[ ! -d $1/RNApaln ]]; then mkdir $1/RNApaln ; fi                 # pre-existing files get removed in the sub-commands
 # if [[ -e ${1}/qsub.RNApaln.out ]] ; then rm ${1}/qsub.RNApaln.out ; fi 
 # if [[ -e ${1}/qsub.RNApaln.err ]] ; then rm ${1}/qsub.RNApaln.err ; fi 
 # CMD="qsub -V -cwd -P ${GRID_ACCOUNT} -N RNApaln -o ${1}/logs/RNApaln.$(date +"%d%m%y-%H%M").out -e ${1}/logs/RNApaln.$(date +"%d%m%y-%H%M").err "
 # # Cannot be parallelized as it reads/writes to flatfiles with static names (1_dp.ps & 2_dp.ps)
 # CMD=${CMD}" -b y ./RNApaln.sge ${1}"
 # RNApaln=$( $CMD ) 
 # echo "launching: "$CMD
 # echo $RNApaln
 # not concurrent ,so just run cleanup in sge script
 # JID=$( echo ${RNApaln} | cut -d " " -f 3 | cut -d "." -f 1 )
 # qsub -hold_jid $JID -V -o ./$1/logs/RNApaln.$(date +"%d%m%y-%H%M").clean -cwd -N clean_RNApaln -pe smp 1 -l mem_requested=1G,h_vmem=1G -j y -b y ./cleanup.sge $1 RNApaln

################################################################
#                         RUN Dynalign_ii
################################################################


#################################################################
#                       post-processing                              


#################################################################
#                             bash
#################################################################

# tar -c $1/*scores.matrix $1/*_time.txt  $1/*structure_ids.txt $1/pairwise_list.txt | gzip > $1_results.tgz

## SCP to local directory? 
# scp $GAMMA://share/ScratchGeneral/marsmi/rfam/*_results.tgz ./
# for f in *.tar; do tar xf $f; done
# for dir in */ ; do for file in $dir*matrix ; do awk '{print $1,$2,$5}' $file > ${file%*matrix}short ; done; done


#################################################################
#                             R
#################################################################

# R CMD BATCH [options] my_script.R [outfile]


