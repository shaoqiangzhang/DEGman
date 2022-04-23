datapreprocess<-function(Express,min.cells=min.cells,log=F){
## preprocess expression data
  #Express=Express[which(rowSums(Express) > 0),]
  if(missing(min.cells)==TRUE){
    Express=Express[apply(Express,1,function(x) sum(x>0)>=3),]
  }else{
    Express=Express[apply(Express,1,function(x) sum(x>0)>=min.cells),]
  }
  if(log==TRUE){
    Express=log2(Express+1)
  }
  return(Express)
}


DEGman<-function(Express,n1,n2){
##Express is the data of expression (genes*cells)
## n1: the cell number of the first group
## n2: the cell number of the second group
  requiredPackages = c('pscl','fitdistrplus','progress')
  for(p in requiredPackages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
  }  
  genelist=rownames(Express)
  Express=apply(Express,2,as.numeric)
  e1=Express[,1:n1]
  e2=Express[,(n1+1):(n1+n2)]
  gene_pval<-c()
  #gene_bhat<-c()
  genename<-c()
  
  pb <- progress_bar$new(
  format = "  DEGman is processing [:bar] :current of total :total genes",
  total = nrow(e1), clear = FALSE, width= 80)

  
  for (i in 1:nrow(e1)) {
    group1=round(as.vector(e1[i,]))
    group2=round(as.vector(e2[i,]))
    bhat=Bhattacharyya_distance(group1,group2)
    
    if(bhat>=0.15){ #default=0.2
      groups=as.vector(Express[i,])
      count=0
      for(times in 1:1000){
      ##random 1000 samples 
      ##find number of samples > original Bhatt distance
        #groups=sample(groups,length(groups))
        groups=sample(groups)
        grp1=groups[1:n1]
        grp2=groups[(n1+1):(n1+n2)]
        if(Bhattacharyya_distance(grp1,grp2)>bhat){
          count=count+1
        }
      }
      pval=count/1000
        if(pval<0.005){#default=0.005
          gene_pval<-c(gene_pval,pval)
          genename<-c(genename,genelist[i])
        }else if(pval<0.1){#default=1
          distribution1=bestFitModel(group1)
          prob1=DistrOfFit(group1,distribution1)
          distribution2=bestFitModel(group2)
          prob2=DistrOfFit(group2,distribution2)
          fitBhat=BhattDistance(prob1,prob2)
          if(fitBhat>=0.15){ #default=0.15
            gene_pval<-c(gene_pval,pval)
            genename<-c(genename,genelist[i])
          }
        }
    }
    pb$tick()
  }

  dataf=data.frame(genename,gene_pval)
  sorted_pval<-dataf[order(dataf$gene_pval),]
  adjusted_pval<-p.adjust(as.numeric(sorted_pval$gene_pval),method="fdr",n=length(as.numeric(sorted_pval$gene_pval)))
  dataf=cbind(sorted_pval,adjusted_pval)
  dataf=dataf[which(dataf$adjusted_pval<0.05),] 
  return(dataf)
}
