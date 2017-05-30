cat("File name","TP","TN","FP","FN","SENS","SPEC","ACC","\n",sep="\t", file="accuracies.tsv")
file.names <- dir(pattern="*_clust.tsv$")
for(x in 1:length(file.names)){
gc <- read.delim(file.names[x], header=F)
# for 1 - max V2
TP=0
FP=0 
NumClust <- max(gc$V2)
for ( cl in 0:NumClust) {
  if ( cl %in% gc$V2 ) {
  v <- as.vector( gc$V1[ gc$V2 == cl ] )
  t <- sort( table( v ), decreasing=T )
  best <- as.integer( t[1] )
  cID <- names( t[ 1 ] )
  if ( cl == 0 ) {
    if ( cID == "shuffled" ) {
      FN <- length(v)-best
      TN <- best
    }
    else 
      cat("Houston, we have a TN problem")
  }
  else {
    if ( cID == "shuffled" ) {
      FP = FP + length(v)
    }
    if ( is.na( as.integer( t[2] )) || as.integer( t[2] ) < best )  {
      TP = TP + best
      FP = FP + length(v)-best
    }
    else if ( as.integer( t[2] ) == best ) {
      # treat both as false positives
      FP = FP + length(v)
    }
  }
}}
TP
TN
FP
FN
SENS=TP / (TP + FN )
SENS
SPEC=TN / ( TN + FP )
SPEC
ACC=(TP + TN) / ( TP + TN + FP + FN )  
ACC
cat(file.names[x],TP,TN,FP,FN,SENS,SPEC,ACC,"\n",sep="\t", file="accuracies.tsv", append=T)
}
