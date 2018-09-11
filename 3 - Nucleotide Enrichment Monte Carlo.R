require(doParallel)

#setwd("") #data directory
dfPreProcess <- function(df, strain='', subs =c(1:8)){
  
  df$series <- ifelse(df$rounds==0, 'No_Realignment', 'Realigned') #assign series
  
  #get sub sums
  df$global <- rowSums(df[,subs])
  
  #prepare data
  df <- df[, c(subs, grep('global', colnames(df)), grep('series', colnames(df)), grep('mini_genome', colnames(df)), grep('Trim_Sequence', colnames(df)), grep('NT_coupure', colnames(df)))] #delimit columns
  #df$mini_genome <- substr(df$mini_genome, 100-10, 100+28) #trims from -10 to 28, which is the range of deviance from TSS allowed by upstream data limits
  df$dangling <- substr(df$Trim_Sequence, nchar(df$Trim_Sequence)-1, nchar(df$Trim_Sequence)-1) #dangling end of primer
  df$end <- substr(df$Trim_Sequence, nchar(df$Trim_Sequence), nchar(df$Trim_Sequence)) #end of primer; directed at 3' U1
  df$upstream <- substr(df$NT_coupure, 1, 1) #1nt upstream of cleavage site
  df <- df[,-c(grep('Trim_Sequence', colnames(df)), grep('NT_coupure', colnames(df)))] #drop useless columns
  df$strain <- strain
  
  gc()
  return(df)
}

names <- c('pr8', 'hk', 'wsn', 'bri')
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
cluster <- ifelse(detectCores()-1 < length(names), detectCores() - 1, length(names))
cl <- makeCluster(cluster)
registerDoParallel(cl)

#parallel import
dfs <- foreach(j = 1:length(names)) %dopar% 
  dfPreProcess(read.csv(paste(toupper(names[j]), '_All_Match.csv', sep=''), stringsAsFactors = F), strains[j])
stopCluster(cl) #end parallel

#assign dfs to strain dfs
for(i in 1:length(names)){
  assign(names[i], dfs[[i]])
}
rm(dfs, strains)
gc()

#you will need to create a p-values folder within the G+1 folder
setwd("./G+1/p-values")
cleavageMC <- function(df, sub=9, iterations=100000, cores=30){
  #takes a data.frame with readcounts as specified in column sub
  # runs a monte carlo procedure on the dangling end, 3' U1 targeted
  # and the nucleotide at +1 of the cleavage site
  #
  #Args:
  # df: data.frame pre-processed for this procedure
  # sub:  column with the subunit of interest
  
  #names for saving .csv
  strain <- df$strain[1]
  subName <- colnames(df)[sub]
  
  #remove rows with 0s
  df <- df[df[,sub]>0,]
  
  #obtain total reads
  allReads <- sum(df[,sub])
  no_RealignmentReads <- sum(df[,sub][df$series=='No_Realignment'])
  realignedReads <- sum(df[,sub][df$series=='Realigned'])
  
  #obtain true values
  #36 numbers:
  #All dangling{A, C, G, U} end{A, C, G, U} upstream{A, C, G, U} | No_Realignment angle{A, C, G, U} end{A, C, G, U} upstream{A, C, G, U} | Realigned angle{A, C, G, U} end{A, C, G, U} upstream{A, C, G, U}
  trueValues <- c(sum(df[,sub][df$dangling=='A'])/allReads * 100,
                sum(df[,sub][df$dangling=='C'])/allReads * 100,
                sum(df[,sub][df$dangling=='G'])/allReads * 100,
                sum(df[,sub][df$dangling=='T'])/allReads * 100,
                sum(df[,sub][df$end=='A'])/allReads * 100,
                sum(df[,sub][df$end=='C'])/allReads * 100,
                sum(df[,sub][df$end=='G'])/allReads * 100,
                sum(df[,sub][df$end=='T'])/allReads * 100,
                sum(df[,sub][df$upstream=='A'])/allReads * 100,
                sum(df[,sub][df$upstream=='C'])/allReads * 100,
                sum(df[,sub][df$upstream=='G'])/allReads * 100,
                sum(df[,sub][df$upstream=='T'])/allReads * 100,
                sum(df[,sub][df$dangling=='A' & df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$dangling=='C'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$dangling=='G'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$dangling=='T'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$end=='A'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$end=='C'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$end=='G'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$end=='T'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$upstream=='A'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$upstream=='C'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$upstream=='G'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$upstream=='T'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                sum(df[,sub][df$dangling=='A' & df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$dangling=='C'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$dangling=='G'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$dangling=='T'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$end=='A' & df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$end=='C'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$end=='G'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$end=='T'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$upstream=='A'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$upstream=='C'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$upstream=='G'& df$series=='Realigned'])/realignedReads * 100,
                sum(df[,sub][df$upstream=='T'& df$series=='Realigned'])/realignedReads * 100)
  
  write.csv(t(trueValues), paste(strain, ' ', subName, ' ', 'obtained', '.csv', sep=''), row.names = F)
  
  #generate data.frames for monto carlo
  df <- df[,c(sub, grep('series', colnames(df)), grep('mini_genome', colnames(df)))]
  gc()
  
  #obtain rows for sampling in monte carlo
  mcNumber <- nrow(df)
  
  #get genome size minus 2
  genomeSize <- nchar(df$mini_genome[1]) - 2
  
  monteCarlo <- function(df, mcNumber, genomeSize, allReads, no_RealignmentReads, realignedReads){
    #obtains nt frequencies at a random point in the passed mini-genome
    #
    #Args:
    # df: data.frame with min_genome, series, and read counts
    # mcNumber: number of rows to pull the nt from
    # genomeSize: length of mini_genome -2
    #
    #return:
    # mcStats:  nt frequences overall and for each series
    
    mcPoint <- sample(genomeSize, mcNumber, replace = T)
    df$dangling <- substr(df$mini_genome, mcPoint, mcPoint)
    df$end <- substr(df$mini_genome, mcPoint + 1, mcPoint + 1)
    df$upstream <- substr(df$mini_genome, mcPoint + 2, mcPoint + 2)
    
    #36 numbers:
    #All dangling{A, C, G, U} end{A, C, G, U} upstream{A, C, G, U} | No_Realignment angle{A, C, G, U} end{A, C, G, U} upstream{A, C, G, U} | Realigned angle{A, C, G, U} end{A, C, G, U} upstream{A, C, G, U}
    mcVector <- c(sum(df[,1][df$dangling=='A'])/allReads * 100,
                  sum(df[,1][df$dangling=='C'])/allReads * 100,
                  sum(df[,1][df$dangling=='G'])/allReads * 100,
                  sum(df[,1][df$dangling=='T'])/allReads * 100,
                  sum(df[,1][df$end=='A'])/allReads * 100,
                  sum(df[,1][df$end=='C'])/allReads * 100,
                  sum(df[,1][df$end=='G'])/allReads * 100,
                  sum(df[,1][df$end=='T'])/allReads * 100,
                  sum(df[,1][df$upstream=='A'])/allReads * 100,
                  sum(df[,1][df$upstream=='C'])/allReads * 100,
                  sum(df[,1][df$upstream=='G'])/allReads * 100,
                  sum(df[,1][df$upstream=='T'])/allReads * 100,
                  sum(df[,1][df$dangling=='A' & df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$dangling=='C'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$dangling=='G'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$dangling=='T'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$end=='A'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$end=='C'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$end=='G'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$end=='T'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$upstream=='A'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$upstream=='C'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$upstream=='G'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$upstream=='T'& df$series=='No_Realignment'])/no_RealignmentReads * 100,
                  sum(df[,1][df$dangling=='A' & df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$dangling=='C'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$dangling=='G'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$dangling=='T'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$end=='A'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$end=='C'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$end=='G'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$end=='T'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$upstream=='A'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$upstream=='C'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$upstream=='G'& df$series=='Realigned'])/realignedReads * 100,
                  sum(df[,1][df$upstream=='T'& df$series=='Realigned'])/realignedReads * 100)
    
    #slows down, but prevents RAM leak
    rm(df, mcNumber, genomeSize, allReads, no_RealignmentReads, realignedReads, mcPoint)
    gc()
    
    return(mcVector)
  }
  
  #commence the parallel...
  #burn the CPU in hellfire
  cl <- makeCluster(cores) 
  registerDoParallel(cl)
  
  #monte carlo procedure data frame
  mcdf <- foreach(j=1:iterations, .combine = rbind, .inorder = F) %dopar%
    monteCarlo(df, mcNumber, genomeSize, allReads, no_RealignmentReads, realignedReads)
  stopCluster(cl)
  
  #name the columns
  colnames(mcdf) <- c("All_Dangling_A", "All_Dangling_C", "All_Dangling_G", "All_Dangling_U", "All_End_A", "All_End_C", "All_End_G", "All_End_U", "All_Upstream_A", "All_Upstream_C", "All_Upstream_G", "All_Upstream_U", "No_Realignment_Dangling_A", "No_Realignment_Dangling_C", "No_Realignment_Dangling_G", "No_Realignment_Dangling_U", "No_Realignment_End_A", "No_Realignment_End_C", "No_Realignment_End_G", "No_Realignment_End_U", "No_Realignment_Upstream_A", "No_Realignment_Upstream_C", "No_Realignment_Upstream_G", "No_Realignment_Upstream_U", "Realigned_Dangling_A", "Realigned_Dangling_C", "Realigned_Dangling_G", "Realigned_Dangling_U", "Realigned_End_A", "Realigned_End_C", "Realigned_End_G", "Realigned_End_U", "Realigned_Upstream_A", "Realigned_Upstream_C", "Realigned_Upstream_G", "Realigned_Upstream_U")
  
  #write the simulated data
  write.csv(mcdf, paste(strain, ' ', subName, ' ', 'simulated', '.csv', sep=''), row.names = F)
  
  #obtain z-scores and p-values
  zScore <- function(trueValue, simulated){
    #converts the trueValue into a z-score
    #
    #args:
    # trueValue  value to obtain a z-score for
    # simulated  simulated data
    #
    #return:
    # z z-statistic
    
    z <- (trueValue - mean(simulated))/sd(simulated) #z-score
    
    gc()
    return(z)
  }
  pScore <- function(trueValue, simulated){
    #converts the trueValue into a p-value
    #
    #args:
    # trueValue  value to obtain a p-value for
    # simulated  simulated data
    #
    #return:
    # p p-value as number of rows above or below the value
    
    simMean <- mean(simulated)
    simIterations <- length(simulated)
    
    if(trueValue >= simMean){
      p <- length(simulated[simulated>=trueValue])/simIterations
    } else{
      p <- length(simulated[simulated<=trueValue])/simIterations
    }
    rm(simMean, simIterations)
    gc()
    return(p)
  }
  
  #setup parallel
  zCores <- ifelse(detectCores()-1 < length(trueValues)/4, detectCores() - 1, length(trueValues)/4)
  cl <- makeCluster(zCores) 
  registerDoParallel(cl)
  
  #obtain z scores
  zVersusApplicable <- foreach(j=1:length(trueValues)) %dopar%
    zScore(trueValues[j], mcdf[,j])
  stopCluster(cl)
  
  zVersusApplicable <- unlist(zVersusApplicable)
  pVersusApplicable <- pnorm(-abs(zVersusApplicable))*2
  
  #pscore
  pCores <- ifelse(detectCores()-1 < length(trueValues)/4, detectCores() - 1, length(trueValues)/4)
  cl <- makeCluster(pCores) 
  registerDoParallel(cl)
  
  #obtain p scores
  pScoreVersusApplicable <- foreach(j=1:length(trueValues)) %dopar%
    pScore(trueValues[j], mcdf[,j])
  stopCluster(cl)
  
  pScoreVersusApplicable <- unlist(pScoreVersusApplicable)
  
  #samething but versus genome
  allCols <- c(1:12, 1:12, 1:12) #force col numbers to always take all data instead of specific
  zCores <- ifelse(detectCores()-1 < length(trueValues)/4, detectCores() - 1, length(trueValues)/4)
  cl <- makeCluster(zCores) 
  registerDoParallel(cl)
  
  zVersusAll <- foreach(j=1:length(trueValues)) %dopar%
    zScore(trueValues[j], mcdf[,allCols[j]])
  stopCluster(cl)
  
  zVersusAll <- unlist(zVersusAll)
  pVersusAll <- pnorm(-abs(zVersusAll))*2
  
  pCores <- ifelse(detectCores()-1 < length(trueValues)/4, detectCores() - 1, length(trueValues)/4)
  cl <- makeCluster(pCores) 
  registerDoParallel(cl)
  
  pScoreVersusAll <- foreach(j=1:length(trueValues)) %dopar%
    pScore(trueValues[j], mcdf[,allCols[j]])
  stopCluster(cl)
  
  pScoreVersusAll <- unlist(pScoreVersusAll)
  
  trueValues <- rbind(
    t(trueValues),
    t(as.vector(colMeans(mcdf))),
    t(as.vector(apply(mcdf, 2, sd))),
    t(zVersusApplicable),
    t(zVersusAll),
    t(pVersusApplicable),
    t(pVersusAll), 
    t(pScoreVersusApplicable),
    t(pScoreVersusAll)
  )
  
  colnames(trueValues) <- c("All_Dangling_A", "All_Dangling_C", "All_Dangling_G", "All_Dangling_U", "All_End_A", "All_End_C", "All_End_G", "All_End_U", "All_Upstream_A", "All_Upstream_C", "All_Upstream_G", "All_Upstream_U", "No_Realignment_Dangling_A", "No_Realignment_Dangling_C", "No_Realignment_Dangling_G", "No_Realignment_Dangling_U", "No_Realignment_End_A", "No_Realignment_End_C", "No_Realignment_End_G", "No_Realignment_End_U", "No_Realignment_Upstream_A", "No_Realignment_Upstream_C", "No_Realignment_Upstream_G", "No_Realignment_Upstream_U", "Realigned_Dangling_A", "Realigned_Dangling_C", "Realigned_Dangling_G", "Realigned_Dangling_U", "Realigned_End_A", "Realigned_End_C", "Realigned_End_G", "Realigned_End_U", "Realigned_Upstream_A", "Realigned_Upstream_C", "Realigned_Upstream_G", "Realigned_Upstream_U")
  rownames(trueValues) <- c('composition', 'simulated_mean', 'simulated_sd', 'zApplicable', 'zAll', 'pApplicable', 'pAll','pScoreApplicable', 'pScoreAll')
  write.csv(trueValues, paste(strain, ' ', subName, ' ', 'obtained z and p', '.csv', sep=''), row.names = T)
}

cleavageMC(wsn)
cleavageMC(pr8)
cleavageMC(hk)
cleavageMC(bri)