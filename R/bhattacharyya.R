Bhattacharyya_distance<-function(groupA,groupB){
## groupA and groupB are two vectors
## calculate the Bhattacharyya distance between groupA and groupB
  groupA=round(groupA)
  groupB=round(groupB)
  maxexp=max(max(groupA),max(groupB))+1  #max expression value +1
  setA<-replicate(maxexp,0)## create new all-zero vectors
  setB<-replicate(maxexp,0)
  for (i in 1: length(groupA)) {
  ##Frequency of expressionz values is counted
    setA[groupA[i]+1]=setA[groupA[i]+1]+1
  }
  setA=setA/sum(setA) 
  for (i in 1:length(groupB)){
    setB[groupB[i]+1]=setB[groupB[i]+1]+1
  }
  setB=setB/sum(setB)
  bhat_coefficient=0
  for (j in 1:maxexp){
    bhat_coefficient=bhat_coefficient+sqrt(setA[j]*setB[j])
  }
  if(bhat_coefficient>1){
    bhat_coefficient=1
  }
  bhat_dist=sqrt(1-bhat_coefficient)
  return(bhat_dist)
}

BhattDistance<-function(prob1,prob2){
##prob1 and prob2 are two vectors of probability distributions

  nc=min(length(prob1),length(prob2))
  bhat_coefficient=0
  for (j in 1:nc){
    bhat_coefficient=bhat_coefficient+sqrt(prob1[j]*prob2[j])
  }
  if(bhat_coefficient>1){
    bhat_coefficient=1
  }
  bhat_dist=sqrt(1-bhat_coefficient)
  return(bhat_dist)
}

#bhatt_distance<-function(Express,n1,n2){
#  genelist=rownames(Express)
#  Express=apply(Express,2,as.numeric)
#  e1=Express[,1:n1]
#  e2=Express[,(n1+1):(n1+n2)]
#  allgenes_bhatt<-c()
#  for (i in 1:nrow(e1)) {
#    group1=as.vector(e1[i,]) 
#    group2=as.vector(e2[i,])
#    bhat=Bhattacharyya_distance(group1,group2)
#    allgenes_bhatt<-c(allgenes_bhatt,bhat)
#  }
#  dataf2=data.frame(genelist,allgenes_bhatt)
#  return(dataf2)
#}
