seqGen <-  function(nts=c('A', 'C', 'G', 'U'), length=5, results=NA){
  #takes a set of nucleotides and returns all possible combinations at a length given
  if(is.na(results[1])){
    results <- nts
  } else {
    tempResults <- results
    results <- vector(mode='character', length=(length(tempResults)*length(nts)))
    k <- 1
    for(i in 1:length(tempResults)){
      for(j in 1:length(nts)){
        results[k] <- paste(tempResults[i], nts[j], sep='')
        k <- k + 1
      }
    }
  }
  
  length <- length - 1
  
  if(length==0){
    return(results)
  } else {
    return(seqGen(nts, length, results))
  }
}

standard <- seqGen(length=3)
gEnd <- seqGen(length=4)

#this generates all possible cleavage sites
#cleavageSites <- c(paste(standard[grep('..[^G]', standard)], '.', sep=''), #no G upstream NN|N.
#                   paste(gEnd[grep('..G.', gEnd)],sep='')) #G to 3'C2; NNG|+N

#.AG. and CAG. are the most important and this delimits the search to them
cleavageSites <- c('.AG.', 'CAG.')

rm(standard, gEnd)

#setwd("") # location of the data
dfPreProcess <- function(df, strain='', subs =c(1:8)){
  #cuts data sets down to subs | all reads | series | mini_genome | cleavage
  
  df$series <- ifelse(df$rounds==0, 'No_Realignment', 'Realigned') #assign series
  
  #get sub sums
  df$global <- rowSums(df[,subs])
  
  #prepare data
  df <- df[, c(subs, grep('global', colnames(df)), grep('series', colnames(df)),
               grep('mini_genome', colnames(df)), grep('Trim_Sequence', colnames(df)),
               grep('NT_coupure', colnames(df)))] #delimit columns
  
  #make T U
  df$Trim_Sequence <- gsub('T', 'U', df$Trim_Sequence)
  df$mini_genome <- gsub('T', 'U', df$mini_genome)
  df$NT_coupure <- gsub('T', 'U', df$NT_coupure)
  
  df$cleavage <- paste(substr(df$Trim_Sequence, nchar(df$Trim_Sequence)-1, nchar(df$Trim_Sequence)),
                       substr(df$NT_coupure, 1,2), sep='')
  
  df <- df[,-c(grep('Trim_Sequence', colnames(df)), grep('NT_coupure', colnames(df)))] #drop useless columns
  df$strain <- strain
  
  gc()
  return(df)
}

require(doParallel)
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

setwd("./G+1/Cleavage Sites") #directory must be created
allCleavageMC <- function(df, cleavageSites, sub=9, iterations=100000, cores=30){
  #takes a data.frame with readcounts as specified in column sub
  # runs a monte carlo procedure on the 3nt cleavage site
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
  
  #trueValues matrix
  trueValues <- matrix(nrow = 8, ncol = length(cleavageSites))
  colnames(trueValues) <- cleavageSites
  
  for(i in 1:ncol(trueValues)){
    trueValues[1,i] <- sum(df[,sub][grepl(colnames(trueValues)[i], df$cleavage)])
    trueValues[2,i] <- sum(df[,sub][grepl(colnames(trueValues)[i], df$cleavage) & df$series=='No_Realignment'])
    trueValues[3,i] <- sum(df[,sub][grepl(colnames(trueValues)[i], df$cleavage) & df$series=='Realigned'])
  }
  trueValues[4,] <- trueValues[1,] / allReads * 100
  trueValues[5,] <- trueValues[2,] / no_RealignmentReads * 100
  trueValues[6,] <- trueValues[3,] / realignedReads * 100
  trueValues[7:8,] <- trueValues[2:3,] / trueValues[c(1,1),] * 100
  
  row.names(trueValues) <- c('All Reads', 'No Realignment Reads', 'Realigned Reads',
                             'All %', 'No Realignment %', 'Realigned %',
                             'No PAR', 'PAR')
  
  write.csv(trueValues, paste(strain, ' ', subName, ' ', 'Cleavage Sites', '.csv', sep=''), row.names = T)
  
  #generate data.frames for monto carlo
  df <- df[,c(sub, grep('series', colnames(df)), grep('mini_genome', colnames(df)))]
  gc()
  
  #obtain rows for sampling in monte carlo
  mcNumber <- nrow(df)
  
  #get genome size minus 4 so things can be selected at least 5 nt long
  genomeSize <- nchar(df$mini_genome[1]) - 3
  
  monteCarlo <- function(df, cleavageSites, mcNumber, genomeSize, allReads){
    #obtains nt frequencies at a random point in the passed mini-genome
    #
    #Args:
    # df: data.frame with min_genome, series, and read counts
    # mcNumber: number of rows to pull the nt from
    # genomeSize: length of mini_genome -3
    #
    #return:
    # mcStats:  nt frequences overall and for each series
    
    mcPoint <- sample(genomeSize, mcNumber, replace = T)
    df$cleavage <- substr(df$mini_genome, mcPoint, mcPoint+3)
    
    mcVector <- vector(mode = 'numeric', length=length(cleavageSites))
    for(i in 1:length(mcVector)){
      mcVector[i] <- sum(df[,1][grepl(cleavageSites[i], df$cleavage)]) / allReads * 100
    }
    
    #slows down, but prevents RAM leak
    rm(df, cleavageSites, mcNumber, genomeSize, allReads, mcPoint)
    gc()
    
    return(mcVector)
  }
  
  #commence the parallel...
  #burn the CPU in hellfire
  cl <- makeCluster(cores) 
  registerDoParallel(cl)
  
  #monte carlo procedure data frame
  mcdf <- foreach(j=1:iterations, .combine = rbind, .inorder = F) %dopar%
    monteCarlo(df, cleavageSites, mcNumber, genomeSize, allReads)
  stopCluster(cl)
  
  #name the columns
  colnames(mcdf) <- cleavageSites
  
  #write the simulated data
  write.csv(mcdf, paste(strain, ' ', subName, ' ', 'simulated', '.csv', sep=''), row.names = F)
  
  #mc statistics
  trueValues <- data.frame(trueValues[4:6,])
  trueValues[nrow(trueValues)+1,] <- colMeans(mcdf)
  row.names(trueValues)[nrow(trueValues)] <- 'Sumulated mean'
  trueValues[nrow(trueValues)+1,] <- apply(mcdf, 2, sd)
  row.names(trueValues)[nrow(trueValues)] <- 'Sumulated sd'
  
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
  
  #z, p, and pMC
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'All z Score'
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'No Realignment z Score'
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'Realigned z Score'
  
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'All p-value'
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'No Realignment p-value'
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'Realigned p-value'
  
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'All p Score'
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'No Realignment p Score'
  trueValues[nrow(trueValues)+1,] <- 0 
  row.names(trueValues)[nrow(trueValues)] <- 'Realigned p Score'
  
  for(i in 1:ncol(trueValues)){
    #zscores
    trueValues[6,i] <- zScore(trueValues[1,i], mcdf[,i])
    trueValues[7,i] <- zScore(trueValues[2,i], mcdf[,i])
    trueValues[8,i] <- zScore(trueValues[3,i], mcdf[,i])
    
    #p-values
    trueValues[9,i] <- pnorm(-abs(trueValues[6,i]))*2
    trueValues[10,i] <- pnorm(-abs(trueValues[7,i]))*2
    trueValues[11,i] <- pnorm(-abs(trueValues[8,i]))*2
    
    #pscores
    trueValues[12,i] <- pScore(trueValues[1,i], mcdf[,i])
    trueValues[13,i] <- pScore(trueValues[2,i], mcdf[,i])
    trueValues[14,i] <- pScore(trueValues[3,i], mcdf[,i])
  }
  
  write.csv(trueValues, paste(strain, ' ', subName, ' ', 'obtained z and p', '.csv', sep=''), row.names = T)
}

allCleavageMC(wsn, cleavageSites)
allCleavageMC(pr8, cleavageSites)
allCleavageMC(hk, cleavageSites)
allCleavageMC(bri, cleavageSites)
