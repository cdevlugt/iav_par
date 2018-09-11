require(doParallel)

subunitDecompress <- function(df, sub=1, drops=c(1:8), drop=T, name="Subunit"){
  # decompresses the data by replicating a column based on read count; can also remove columns
  #
  # Args:
  #   df: data frame containing data
  #   sub:  column containing data to be replicated
  #   drops:  columns to be dropped
  #   drop: should the other subunits be dropped
  #
  # Returns:
  #   df: decompress data frame for a given subunit
  if(drop){
    drops <- drops[drops != sub]#retains subunit column if in drops range
    if(length(drops)!=0){
      df <- df[,-drops]
      sub <- sub - sum(drops < sub)#one more opperation versus replicating than dropping, but stores less stuff in ram during replication
    }
  }
  df <- df[(rep(row.names(df), df[,sub])),]
  df[,sub] <- colnames(df)[sub]
  colnames(df)[sub] <- name
  return (df)
}
subunitSeperate <- function(df, subs=c(1:8), decompress=F, name="Subunit"){
  #Takes a data frame with read counts in columns and then clones the rows seperating the subunits so that they are on seperate rows
  #
  #Args:
  # df: data frame with subunit read counts
  # subs: column numbers of subunits
  # decompress: flag indicating if the subunits should be decompressed
  # name: name for collumn if decompress is T
  #
  #Return:
  # dfFinal: data frame with columns cloned based on subunit
  
  for(i in 1:length(subs)){
    dfTemp <- df[df[,subs[i]]>0,]
    if(decompress){
      dfTemp <- subunitDecompress(dfTemp, sub=subs[i], drops=subs, drop=T, name=name)
    } else {
      dfTemp[,subs[-i]] <- 0
    }
    if(i==1){
      dfFinal <- dfTemp
    } else {
      dfFinal <- rbind(dfFinal, dfTemp)
    }
    rm(dfTemp)
  }
  return(dfFinal)
}
parAdditionFreq <- function(df, strain = ''){
  #subunit handling from strain... lets me doParallel the pre-processing
  if(strain=='Luciferase') {
    U = c(2,4)
    C = c(1,3)
  } else if (strain=='Hong Kong') {
    U = c(1:6)
    C = c(7:8)
  } else {
    U = c(1:5)
    C = c(6:8)
  }
  
  U4 <- colnames(df)[U]
  C4 <- colnames(df)[C]
  
  #obtain data needed
  df <- df[,c(U, C, grep('Seq_Trim_R', colnames(df)))]
  gc()
  
  #drop rows with no PAR
  df <- df[df$Seq_Trim_R1!='',]
  
  #expand data
  df <- subunitSeperate(df, c(U, C), T)
  
  #PAR additions
  additions <- c('G', 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA')
  
  subs <- unique(df$Subunit)
  rounds <- length(grep('Seq_Trim_R', colnames(df)))
  
  additionsTable <- function(df, sub='All', rounds=1, additions=c('G', 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA')){
    
    #generate return table
    adt <- data.frame(target=sub, additions=additions, stringsAsFactors = F)
    
    #changing all values is faster than two if statements
    if(sub=='All'){
      df$Subunit <- 'All'
    }
    print(sub)
    df <- df[df$Subunit==sub,]
    gc()
    
    for(i in 1: rounds){
      #create a new column
      adt[, (ncol(adt)+1)] <- 0
      colnames(adt)[ncol(adt)] <- paste('Round_', i, sep='')
      
      for(j in 1: length(additions)){ #populate column
        adt[j, ncol(adt)] <- nrow(df[df[,(i+1)]==additions[j],])
      }
    }
    
    #total for each addition and percent of all PAR
    adt$total <- rowSums(adt[,3:ncol(adt)])
    adt$percent <- adt$total/sum(adt$total) * 100
    
    return(adt)
  }
  
  for(i in 0:length(subs)){
    if(i==0){
      adt <- additionsTable(df, 'All', rounds, additions)
    } else {
      adt <- rbind(adt, additionsTable(df, subs[i], rounds, additions))
    }
  }
  
  #templates handling
  df$Subunit[df$Subunit %in% U4] <- "3'U4"
  df$Subunit[df$Subunit %in% C4] <- "3'C4"
  
  subs <- c("3'U4", "3'C4")
  
  for(i in 1:length(subs)){
    adt <- rbind(adt, additionsTable(df, subs[i], rounds, additions))
  }
  setwd("./Sequence Additions") #directory must exist
  write.csv(adt, paste(strain, ' nucleotide sequence additions.csv', sep=''))
}

#setwd("")
names <- c('pr8', 'hk', 'wsn', 'bri')
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
cluster <- length(names)
cl <- makeCluster(cluster)
registerDoParallel(cl)

dfs <- foreach(j = 1:cluster, .inorder = F) %dopar% 
  parAdditionFreq(read.csv(paste(toupper(names[j]), '_All_Match.csv', sep=''), stringsAsFactors = F), strains[j])
stopCluster(cl)
gc()