#requires 5 - Sequence Addition Table
#####Correlations between length of the Primer and length of the Addition
require(compiler)
enableJIT(3)

require(ggplot2)
require(ggpubr)
require(cowplot)

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
    df <- df[,-drops]
    sub <- sub - sum(drops < sub)#one more opperation versus replicating than dropping, but stores less stuff in ram during replication
  }
  df <- df[(rep(row.names(df), df[,sub])),]
  df[,sub] <- colnames(df)[sub]
  colnames(df)[sub] <- name
  return (df)
}
adjustByLength <- function(df, sub=NA, startLen=26, additions=c(29:36), rounds=NA, minAdjust=1, end=F, decompress=F, byAddition=T, byLength=T){
  #takes a data frame containing a start point and lengths of adjustments, and returns a data frame containing
  # the average adjustment either to or from that length. both strings and intger additions may be passed
  #
  #Args:
  # df: the dataframe
  # sub:  the column containing read counts; if unspecified it will take all data
  # startLen: column containing start length
  # additions:  columns containing additions
  # rounds: column containing rounds of addition
  # minAdjust:  minimum adjustment for tabulated data
  # end:  should data be returned as end lengths and average addition to that end length (F)
  #   or as the start length plus average addition from that length(T)
  # decompress: decompresses based on reads in subunit... this requires a large amount of RAM
  # byAddition:  average uncompressed data by addition; if both byAddition and byLength, averages are calculated by global
  # byLength: average uncompressed data by length; if both byAddition and byLength, averages are calculated by global
  #
  #return:
  # dfAdj:  data frame containing lengths and average additions to or from that point depending on specification
  # dfABL:  decompressed file
  
  #determines addition type
  if(is.numeric(df[1,additions[1]])){
    char <- F
  } else if(is.character(df[1,additions[1]])) {
    char <- T
    seqs <- unique(df[,additions[1]])
    if(length(additions) > 1){
      for(i in 2:length(additions)){
        seqs <- c(seqs, df[,additions[i]])
      }
    }
    seqs <- unique(seqs)
    seqs <- seqs[seqs!=""]
    seqs <- seqs[order(nchar(seqs), seqs)]
  } else {
    stop("Invalid type in additions argument")
  }
  
  #if specified removes all rows with rounds of 0
  if(!(is.na(rounds))){
    df <- df[df[,rounds]>0,]
  } else {
    df[,(ncol(df)+1)] <- 0
    rounds <- ncol(df)
  }
  
  #removes rows with less than the minimum adjustment
  if(!char){
    df <- df[rowSums(df[,additions])>=minAdjust,]
  } else {
    df$tempcol <- df[,additions[1]]
    if(length(additions)>1){
      for(i in 2:length(additions)){
        df$tempcol <- paste(df$tempcol, df[,additions[i]], sep="")
      }
    }
    df <- df[nchar(df$tempcol)>=minAdjust,]
    df <- df[,-ncol(df)]
  }
  
  #takes all if sub unspecified, drops all with sub reads == 0
  if(is.na(sub)){
    df[,(ncol(df)+1)] <- 1
    sub <- ncol(df)
  }
  df <- df[df[,sub]>0,]
  
  #removes useless data
  cols <- c(sub, rounds, startLen, additions)
  df <- df[,cols]
  gc()
  
  #gets length adjustments 
  dfAdj <- df[, c(1, 3, 4)]
  if(!char){
    dfAdj[,(ncol(dfAdj)+1)] <- dfAdj[,2] + dfAdj[,3]
  } else {
    dfAdj[,(ncol(dfAdj)+1)] <- dfAdj[,2] + nchar(dfAdj[,3])
  }
  
  #obtain remainder of additions
  if(length(additions)>2){
    for(i in 2:length(additions)){
      dfAdj[,(ncol(dfAdj)+1)] <- df[,(3+i)]
      if(!char){
        dfAdj[,(ncol(dfAdj)+1)] <- dfAdj[,ncol(dfAdj)] + dfAdj[,(ncol(dfAdj)-1)]
        dfAdj[,ncol(dfAdj)][dfAdj[,(ncol(dfAdj)-1)]==0] <- 0
      } else {
        dfAdj[,(ncol(dfAdj)+1)] <- nchar(dfAdj[,ncol(dfAdj)]) + dfAdj[,(ncol(dfAdj)-1)]
        dfAdj[,ncol(dfAdj)][dfAdj[,(ncol(dfAdj)-1)]==""] <- 0
      }
    }
  }
  rm(df) #removes df
  gc()
  
  #creates data frame as: length, addition, end
  dfABL <- dfAdj[,c(1,2,3)]
  colnames(dfABL)[2:3] <- c("Length", "Addition")
  if(!char){
    dfABL$End <- dfABL$Length + dfABL$Addition
  } else {
    dfABL$End <- dfABL$Length + nchar(dfABL$Addition)
  }
  
  #add in remainder of rounds via rbind
  if(length(additions)>1){
    for(i in 2:length(additions)){
      dfABLtemp <- dfAdj[,c(1,(i*2),(i*2+1))]
      colnames(dfABLtemp)[2:3] <- c("Length", "Addition")
      if(!char){
        dfABLtemp$End <- dfABLtemp$Length + dfABLtemp$Addition
        dfABLtemp <- dfABLtemp[dfABLtemp$Addition > 0,]
      } else {
        dfABLtemp$End <- dfABLtemp$Length + nchar(dfABLtemp$Addition)
        dfABLtemp <- dfABLtemp[nchar(dfABLtemp$Addition) > 0,]
      }
      dfABL <- rbind(dfABL, dfABLtemp)
    }
  }
  rm(dfABLtemp)
  gc()
  
  if(!decompress){ #average by multiplying by read count
    if(!char){
      dfABL[,3] <- dfABL[,3] * dfABL[,1]
      
      if(end){
        dfAdj <- data.frame(table(dfABL$End))
        col <- 4
      } else {
        dfAdj <- data.frame(table(dfABL$Length))
        col <- 2
      }
      dfAdj[,1] <- as.numeric(as.character(dfAdj[,1]))
      
      for(i in 1:nrow(dfAdj)){
        dftemp <- dfABL[dfABL[,col]==dfAdj[i,1],]
        dfAdj[i,2] <- sum(dftemp[,3])/sum(dftemp[,1]) 
      }
    } else {
      if(end){
        col <- 4
      } else {
        col <- 2
      }
      
      #lengths to create data frame
      #lengths <- unique(dfABL[,col])
      #lengths <- lengths[order(lengths)]
      lengths <- 9:16 #data range for all analysis
      dfABL <- dfABL[dfABL$Length>=min(lengths) & dfABL$Length<=max(lengths),]
      
      dfAdj <- data.frame(Length=rep(lengths, length(seqs)), stringsAsFactors = F) #reassigns dfAdj for return
      dfAdj$Addition <- ""
      dfAdj <- dfAdj[order(dfAdj[,1]),]
      dfAdj$Frequency <- 0
      dfAdj$Seq_Reads <- 0
      dfAdj$Reads <- 0
      row.names(dfAdj) <- 1:nrow(dfAdj)
      
      globalCount <- sum(dfABL[,1])
      
      if(byAddition == byLength){
        count <- globalCount
      }
      
      dfAdj$LenCount <- 0
      dfAdj$AddCount <- 0
      dfAdj$GlobalCount <- globalCount
      
      rm(globalCount)
      
      for(i in 1:length(lengths)){
        dfi <- dfABL[dfABL[,col]==lengths[i],]
        lengthCount <- sum(dfi[,1])
        if(!byAddition && byLength){
          count <- lengthCount
        }
        for(j in 1:length(seqs)){
          dfj <- dfi[dfi[,3]==seqs[j],] #data frame with specific length (start or end) for a given addition
          addCount <- sum(dfABL[,1][dfABL[,3]==seqs[j]])
          if(byAddition && !byLength){
            count <- addCount
          }
          dfAdj[((i-1)*length(seqs)+j),2] <- seqs[j] #assigns sequence to Addition column
          dfAdj[((i-1)*length(seqs)+j),3] <- sum(dfj[,1])/count #assigns frequency to frequency column
          dfAdj[((i-1)*length(seqs)+j),4] <- sum(dfj[,1]) #total number of reads
          dfAdj[((i-1)*length(seqs)+j),5] <- count #assigns number of reads of that length to Reads
          dfAdj[((i-1)*length(seqs)+j),6] <- lengthCount
          dfAdj[((i-1)*length(seqs)+j),7] <- addCount
          rm(dfj) #close to save RAM
          gc()
        }
        rm(dfi) #close to save RAM
        gc()
      }
    }
    dfAdj[is.na(dfAdj)] <- 0
    
    return(dfAdj) #heatmap return
  } else { #decompressed table
    if(end){
      col <- 4
    } else {
      col <- 2
    }
    rm(dfAdj)
    gc()
    
    dfABL <- dfABL[,c(1,col,3)] #get reads, start or end, adjustment
    dfABL <- dfABL[order(dfABL[,2], nchar(dfABL[,3]), dfABL[,3]),]
    
    #length delimiting
    dfABL <- dfABL[dfABL[,2] >= 8 & dfABL[,2] <= 17,]
    
    dfABL <- subunitDecompress(dfABL, sub=1, drop=F)
    dfABL <- dfABL[,c(2,3)]
    colnames(dfABL) <- c("length", "addition")
    
    gc()
    return(dfABL)
  }
}
#setwd("") #set directory here
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

setwd("./Sequence Additions") #folder must exist
maxRounds <- max(hk$rounds, pr8$rounds, bri$rounds, wsn$rounds)

#create a large data frame for all strains
pr8$PuertoRico8 <- rowSums(pr8[,c(1:8)])
pr8$HongKong <- 0
pr8$WSN <- 0
pr8$Brisbane <- 0
pr8 <- pr8[,c((ncol(pr8)-3):(ncol(pr8)), grep("Trim_Len", colnames(pr8)), grep("rounds", colnames(pr8)), grep("Num_Trim", colnames(pr8)))]
dfRounds <- max(pr8$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    pr8[,(ncol(pr8)+1)] <- 0
    colnames(pr8)[ncol(pr8)] <- paste("Num_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

hk$PuertoRico8 <- 0
hk$HongKong <- rowSums(hk[,c(1:8)])
hk$WSN <- 0
hk$Brisbane <- 0
hk <- hk[,c((ncol(hk)-3):(ncol(hk)), grep("Trim_Len", colnames(hk)), grep("rounds", colnames(hk)), grep("Num_Trim", colnames(hk)))]
dfRounds <- max(hk$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    hk[,(ncol(hk)+1)] <- 0
    colnames(hk)[ncol(hk)] <- paste("Num_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

wsn$PuertoRico8 <- 0
wsn$HongKong <- 0
wsn$WSN <- rowSums(wsn[,c(1:8)])
wsn$Brisbane <- 0
wsn <- wsn[,c((ncol(wsn)-3):(ncol(wsn)), grep("Trim_Len", colnames(wsn)), grep("rounds", colnames(wsn)), grep("Num_Trim", colnames(wsn)))]
dfRounds <- max(wsn$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    wsn[,(ncol(wsn)+1)] <- 0
    colnames(wsn)[ncol(wsn)] <- paste("Num_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

bri$PuertoRico8 <- 0
bri$HongKong <- 0
bri$WSN <- 0
bri$Brisbane <- rowSums(bri[,c(1:8)])
bri <- bri[,c((ncol(bri)-3):(ncol(bri)), grep("Trim_Len", colnames(bri)), grep("rounds", colnames(bri)), grep("Num_Trim", colnames(bri)))]
dfRounds <- max(bri$rounds)
if(dfRounds < maxRounds){
  for(i in (dfRounds+1):maxRounds){
    bri[,(ncol(bri)+1)] <- 0
    colnames(bri)[ncol(bri)] <- paste("Num_Trim_R", i, sep="")
  }
}
rm(dfRounds)
gc()

df <- rbind(pr8, hk, wsn, bri)
rm(pr8, bri, hk, wsn)
gc()

require(doParallel)
cores <- ifelse(detectCores()<8, detectCores(), 8)
cl <- makeCluster(cores) 
registerDoParallel(cl)

#dopar args
subs <- c(1,1,2,2,3,3,4,4) #each sub twice
ends <- c(F,T,F,T,F,T,F,T) #ends argument

dfs.adj <- foreach(j = 1:8) %dopar%
  adjustByLength(df, sub=subs[j], startLen=grep('Trim_Len', colnames(df)), additions=grep('Trim_R', colnames(df)), rounds=grep('rounds', colnames(df)), minAdjust=1, end=ends[j], decompress=T, byAddition=T, byLength=T)

dfs.pearson <- foreach(j = 1:8) %dopar%
  cor.test(dfs.adj[[j]]$length, dfs.adj[[j]]$addition, method='pearson')

dfs.spearman <- foreach(j = 1:8) %dopar%
  cor.test(dfs.adj[[j]]$length, dfs.adj[[j]]$addition, method='spearman')
stopCluster(cl)


for(i in 1:8){
  dfs.adj[[i]] <- dfs.adj[[i]][dfs.adj[[i]]$addition>=4,]
}

spear.av <- mean(c(dfs.spearman[[1]]$estimate, dfs.spearman[[3]]$estimate, dfs.spearman[[5]]$estimate, dfs.spearman[[7]]$estimate))
spear.sd <- sd(c(dfs.spearman[[1]]$estimate, dfs.spearman[[3]]$estimate, dfs.spearman[[5]]$estimate, dfs.spearman[[7]]$estimate))

pear.av <- mean(c(dfs.pearson[[1]]$estimate, dfs.pearson[[3]]$estimate, dfs.pearson[[5]]$estimate, dfs.pearson[[7]]$estimate))
pear.sd <- sd(c(dfs.pearson[[1]]$estimate, dfs.pearson[[3]]$estimate, dfs.pearson[[5]]$estimate, dfs.pearson[[7]]$estimate))

cdf <- data.frame(name=c(colnames(df[,1:4]), 'mean', 'sd'), 
                  spearman=c(dfs.spearman[[1]]$estimate, dfs.spearman[[3]]$estimate, dfs.spearman[[5]]$estimate, dfs.spearman[[7]]$estimate, spear.av, spear.sd),
                  pearson=c(dfs.pearson[[1]]$estimate, dfs.pearson[[3]]$estimate, dfs.pearson[[5]]$estimate, dfs.pearson[[7]]$estimate, pear.av, pear.sd))

write.csv(cdf, 'Sequence Addition and Length Correlations.csv', row.names = F)

#####Correlations between length of the addition and frequency of the addition
require(doParallel)
require(ggplot2)
require(cowplot)

dfPreProcess <- function(strain = '', cols = c('additions', 'total', 'percent'), type = 'All'){
  #processes the data to length of addition, total reads, percent, ln(percent)
  # strain: reads data from strain name
  # cols: vector of column names
  # type: subunit or template or all to anaylze
  
  #read data
  df <- read.csv(paste(strain, ' nucleotide sequence additions', '.csv', sep=''), stringsAsFactors = F, row.names = 1)
  df <- df[df$target==type,] #limits data to applicable
  df <- df[,grep(paste(cols, collapse = '|'), colnames(df))] #remove useless data
  
  #convert data to lengths rather than sequence
  df.temp <- data.frame(length=1:6, total = 0, percent = 0)
  for(i in 1:6){
    df.temp[i,2] <- sum(df[,2][nchar(df$additions)==i])
  }
  df.temp[,3] <- df.temp[,2]/sum(df.temp[,2])*100
  df.temp$log <- log(df.temp[,3]) #get that log
  df.temp$strain <- strain
  
  return(df.temp)
}

#read data in parallel
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
cluster <- 4
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster, .combine = 'rbind') %dopar%
  dfPreProcess(strains[j])
stopCluster(cl)
gc()

#factor strain for ggplot colour consistency
dfs$strain <- factor(dfs$strain, levels = strains)

#pass data for 3:6 to all for correlations
all <- dfs[dfs$length>2,]

r2val <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  r2val[i] <- summary(lm(data = all[all$strain==strains[i],], log ~ length))$r.squared
}

corPlot <- ggplot(dfs, aes(x = length, y = log, group=strain, colour=strain)) +
  geom_point(size=0.5) +
  geom_smooth(data=all, aes(x = length, y = log, group=strain, colour=strain), method='lm', se=F, formula = y~x, size=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=4),
        axis.text = element_text(size=4), 
        axis.title = element_text(size=6))+
  scale_color_manual(breaks=strains,
                     values=c(hcl(h=seq(15,375, length=(4+1))[1], c=100, l=65), #retain colour
                              hcl(h=seq(15,375, length=(4+1))[2], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[3], c=100, l=65),
                              hcl(h=seq(15,375, length=(4+1))[4], c=100, l=65)),
                     labels=c(do.call('expression', #add R^2 values with supersctipts to legend
                                      list(
                                        bquote(.(strains[1])~(R^{2}*~'='~.(formatC(r2val[1], digits = 3, flag='#')))),
                                        bquote(.(strains[2])~(R^{2}*~'='~.(formatC(r2val[2], digits = 3, flag='#')))),
                                        bquote(.(strains[3])~(R^{2}*~'='~.(formatC(r2val[3], digits = 3, flag='#')))),
                                        bquote(.(strains[4])~(R^{2}*~'='~.(formatC(r2val[4], digits = 3, flag='#'))))
                                      )))) +
  labs(x='Sequence Addition Length (nt)', y='Addition Frequency (ln(%))') +
  scale_x_continuous(limits = c(1, 6), breaks=1:6, labels = 1:6)

#steal legend
corLegend <- get_legend(corPlot)
#remove legend
corPlot <- corPlot + theme(legend.position = 'none')

#save plot
save_plot('Sequence Addition and ln of Frequency.png', plot=plot_grid(corPlot, corLegend, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)

#save r2
write.csv(data.frame(Strain = strains, R2 = r2val, stringsAsFactors = F), file='Sequence addition length and Frequency between 3 and 6 nt.csv')
