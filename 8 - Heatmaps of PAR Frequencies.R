#requires: 8 - Effect Modeling Probability Tables
require(ggplot2)
require(cowplot)
require(doParallel)

#import strain PAR
#setwd() #if not default directory, set it here
setwd('./Overview')
parav <- read.csv('PAR rates.csv', stringsAsFactors = F, row.names = 1)[2,]
stPAR <- as.numeric(rowMeans(parav[1,c(1,2,4)])) #because WSN length is invalid it is omitted

#splits for graph
splits <- c(0, mean(c(mean(c(stPAR, 0)),0)), mean(c(mean(c(stPAR, 0)),stPAR)), #below
            stPAR, #PAR rate
            mean(c(mean(c(stPAR, 100)),stPAR)), mean(c(mean(c(stPAR, 100)),100)), 100) #above
splits <- splits/100 #fix for heatmap

#colour scheme
#this is set up so that the colours are similar to the distance from the mean PAR rate of all 3 strains 
colours = c('green', '#40E0D0', 'light grey',
            'white',
            'yellow', 'orange', 'red')

setwd('..')
setwd('./Effect Modeling')

#####All data heatmaps
#read data in parallel
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
series <- c('Realigned', 'No_Realignment')
cluster <- 8
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[floor((j-1)/2)+1], ' ', series[floor(j/2)-ceiling(j/2)+2], '_counts.csv', sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

#unlist data
varNames <- c('pre', 'pno', 'pal',
              'hre', 'hno', 'hal',
              'wre', 'wno', 'wal',
              'bre', 'bno', 'bal')
for(j in 1:4){
  assign(varNames[(j-1)*3+1], dfs[[(j-1)*2+1]])
  assign(varNames[(j-1)*3+2], dfs[[(j-1)*2+2]])
  assign(varNames[(j-1)*3+3],(dfs[[(j-1)*2+1]] + dfs[[(j-1)*2+2]]))
}

#threshold
#6400th of reads as threshold
threshold <- 6400
pal[pal<=sum(pal)/threshold] <- 0
hal[hal<=sum(hal)/threshold] <- 0
wal[wal<=sum(wal)/threshold] <- 0
bal[bal<=sum(bal)/threshold] <- 0

#obtain percentages
pr8 <- as.matrix(pre/pal*100)
hk <- as.matrix(hre/hal*100)
wsn <- as.matrix(wre/wal*100)
bri <- as.matrix(bre/bal*100)

#drop inf
for(i in 1:ncol(pr8)){
  pr8[,i][is.infinite(pr8[,i])] <- NaN
  hk[,i][is.infinite(hk[,i])] <- NaN
  wsn[,i][is.infinite(wsn[,i])] <- NaN
  bri[,i][is.infinite(bri[,i])] <- NaN
}

#create a matrix with the average of all matrix in varNames
#sums ignoring NaN
varNames <- c('pr8', 'hk', 'bri')
all <- matrix(NaN, nrow = nrow(pr8), ncol = ncol(pr8))
row.names(all) <- row.names(pr8)
colnames(all) <- colnames(pr8)
for(i in 1:ncol(pr8)){
  for(j in 1:nrow(pr8)){
    val <- 0 #placeholder value
    cnt <- 0 #count of vals
    for(k in 1:length(varNames)){
      if(!is.na(get(varNames[k])[j,i])){
        val <- val + get(varNames[k])[j,i]
        cnt <- cnt + 1
      }
    }
    if(cnt>0){
      all[j,i] <- val/cnt
    }
  }
}

#make data ggplot compatible
Seq_Ends <- rep(row.names(all), ncol(all))
lengths <- rep(colnames(all)[1], nrow(all))
PAR <- all[,1]
for(i in 2:ncol(all)){
  lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
  PAR <- c(PAR, all[,i])
}

#put ggplot friendly data into a df
hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
rm(PAR, Seq_Ends, lengths)
gc()

#heatmap
hmA <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
  geom_tile(colour="white",size=0.25) +
  geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
  scale_fill_gradientn(limits=c(0,100),
                       colours = colours,
                       values = splits,
                       na.value='#000000') +
  theme_bw() +
  ggtitle("All Conserved Sequences") +
  theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
  scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
  xlab('Length (nt)') +
  ylab("Sequence End") +
  theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
  theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) +
  theme(legend.title =element_blank(), legend.text=element_text(size=7))
legend <- get_legend(hmA)
hmA <- hmA + theme(legend.position='none')

#####3'U4 data heatmaps
#read data in parallel
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
series <- c('Realigned', 'No_Realignment')
cluster <- 8
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[floor((j-1)/2)+1], ' ', series[floor(j/2)-ceiling(j/2)+2], "_counts_3'U4.csv", sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

#unlist data
varNames <- c('pre', 'pno', 'pal',
              'hre', 'hno', 'hal',
              'wre', 'wno', 'wal',
              'bre', 'bno', 'bal')
for(j in 1:4){
  assign(varNames[(j-1)*3+1], dfs[[(j-1)*2+1]])
  assign(varNames[(j-1)*3+2], dfs[[(j-1)*2+2]])
  assign(varNames[(j-1)*3+3],(dfs[[(j-1)*2+1]] + dfs[[(j-1)*2+2]]))
}

#threshold
#1200th of reads as threshold
threshold <- 1200
pal[pal<=sum(pal)/threshold] <- 0
hal[hal<=sum(hal)/threshold] <- 0
wal[wal<=sum(wal)/threshold] <- 0
bal[bal<=sum(bal)/threshold] <- 0

#calculate PAR rates
pr8 <- as.matrix(pre/pal*100)
hk <- as.matrix(hre/hal*100)
wsn <- as.matrix(wre/wal*100)
bri <- as.matrix(bre/bal*100)

#drop inf
for(i in 1:ncol(pr8)){
  pr8[,i][is.infinite(pr8[,i])] <- NaN
  hk[,i][is.infinite(hk[,i])] <- NaN
  wsn[,i][is.infinite(wsn[,i])] <- NaN
  bri[,i][is.infinite(bri[,i])] <- NaN
}

#create a matrix with the average of all matrix in varNames
#sums ignoring NaN
varNames <- c('pr8', 'hk', 'bri')
all <- matrix(NaN, nrow = nrow(pr8), ncol = ncol(pr8))
row.names(all) <- row.names(pr8)
colnames(all) <- colnames(pr8)
for(i in 1:ncol(pr8)){
  for(j in 1:nrow(pr8)){
    val <- 0 #placeholder value
    cnt <- 0 #count of vals
    for(k in 1:length(varNames)){
      if(!is.na(get(varNames[k])[j,i])){
        val <- val + get(varNames[k])[j,i]
        cnt <- cnt + 1
      }
    }
    if(cnt>0){
      all[j,i] <- val/cnt
    }
  }
}

#make data ggplot compatible
Seq_Ends <- rep(row.names(all), ncol(all))
lengths <- rep(colnames(all)[1], nrow(all))
PAR <- all[,1]
for(i in 2:ncol(all)){
  lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
  PAR <- c(PAR, all[,i])
}

#tabulate ggplot data
hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
rm(PAR, Seq_Ends, lengths)
gc()

hmU <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
  geom_tile(colour="white",size=0.25) +
  geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
  scale_fill_gradientn(colours = colours,
                       values = splits,
                       na.value='#000000') +
  theme_bw() +
  ggtitle("3'U4 Conserved Sequences") +
  theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
  scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
  xlab('Length (nt)') +
  ylab("Sequence End") +
  theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
  theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) + 
  theme(legend.title =element_blank(), legend.text=element_text(size=7))
hmU <- hmU + theme(legend.position='none')

#####3'C4 data heatmaps
#read data in parallel
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
series <- c('Realigned', 'No_Realignment')
cluster <- 8
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[floor((j-1)/2)+1], ' ', series[floor(j/2)-ceiling(j/2)+2], "_counts_3'C4.csv", sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

#unlist data
varNames <- c('pre', 'pno', 'pal',
              'hre', 'hno', 'hal',
              'wre', 'wno', 'wal',
              'bre', 'bno', 'bal')
for(j in 1:4){
  assign(varNames[(j-1)*3+1], dfs[[(j-1)*2+1]])
  assign(varNames[(j-1)*3+2], dfs[[(j-1)*2+2]])
  assign(varNames[(j-1)*3+3],(dfs[[(j-1)*2+1]] + dfs[[(j-1)*2+2]]))
}

#threshold
#12800th of reads as threshold
threshold <- 12800
pal[pal<=sum(pal)/threshold] <- 0
hal[hal<=sum(hal)/threshold] <- 0
wal[wal<=sum(wal)/threshold] <- 0
bal[bal<=sum(bal)/threshold] <- 0

#calculate PAR frequencies
pr8 <- as.matrix(pre/pal*100)
hk <- as.matrix(hre/hal*100)
wsn <- as.matrix(wre/wal*100)
bri <- as.matrix(bre/bal*100)

#drop inf
for(i in 1:ncol(pr8)){
  pr8[,i][is.infinite(pr8[,i])] <- NaN
  hk[,i][is.infinite(hk[,i])] <- NaN
  wsn[,i][is.infinite(wsn[,i])] <- NaN
  bri[,i][is.infinite(bri[,i])] <- NaN
}

#create a matrix with the average of all matrix in varNames
#sums ignoring NaN
varNames <- c('pr8', 'hk', 'bri')
all <- matrix(NaN, nrow = nrow(pr8), ncol = ncol(pr8))
row.names(all) <- row.names(pr8)
colnames(all) <- colnames(pr8)
for(i in 1:ncol(pr8)){
  for(j in 1:nrow(pr8)){
    val <- 0 #placeholder value
    cnt <- 0 #count of vals
    for(k in 1:length(varNames)){
      if(!is.na(get(varNames[k])[j,i])){
        val <- val + get(varNames[k])[j,i]
        cnt <- cnt + 1
      }
    }
    if(cnt>0){
      all[j,i] <- val/cnt
    }
  }
}

#make data ggplot compatible
Seq_Ends <- rep(row.names(all), ncol(all))
lengths <- rep(colnames(all)[1], nrow(all))
PAR <- all[,1]
for(i in 2:ncol(all)){
  lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
  PAR <- c(PAR, all[,i])
}

#tabulate
hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
rm(PAR, Seq_Ends, lengths)
gc()

hmC <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
  geom_tile(colour="white",size=0.25) +
  geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
  scale_fill_gradientn(colours = colours,
                       values = splits,
                       na.value='#000000') +
  theme_bw() +
  ggtitle("3'C4 Conserved Sequences") +
  theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
  scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
  xlab('Length (nt)') +
  ylab("Sequence End") +
  theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
  theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) + 
  theme(legend.title =element_blank(), legend.text=element_text(size=7))
hmC <- hmC + theme(legend.position='none')

#####Save Plots
hm <- plot_grid(hmA, hmU, hmC, legend, ncol=4, nrow=1, rel_widths = c(1,1,1,0.3))
save_plot('Effect Modeling Heatmap.png', hm, base_height = 2.3, base_width = 8, dpi=900)

