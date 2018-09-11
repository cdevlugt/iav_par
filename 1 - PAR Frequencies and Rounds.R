require(doParallel)
require(ggplot2)
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
dfPreProcess <- function(df, strain='', U=1:5, C=6:8) {
  
  #subunit handling from strain... lets me doParallel the pre-processing
  if(is.na(U) && is.na(C) && strain=='Luciferase') {
    U = c(2,4)
    C = c(1,3)
  } else if (is.na(U) && is.na(C) && strain=='Hong Kong') {
    U = c(1:6)
    C = c(7:8)
  } else if (is.na(U) && is.na(C)) {
    U = c(1:5)
    C = c(6:8)
  }
  
  #delimit the data
  df <- df[,c(U, C, grep('rounds', colnames(df)), grep('Seq_Trim_R', colnames(df)))]
  
  #series handing
  df$series <- 'No_Realignment'
  df$series[df$rounds>=1] <- 'Realigned'
  df$series <- factor(df$series, levels=c('Realigned', 'No_Realignment'))
  
  #subunit handling
  subs <- c(U, C)
  subs <- unique(subs)
  subs <- subs[order(subs)]
  subNames <- colnames(df)[subs]
  C <- colnames(df)[C]
  
  df <- subunitSeperate(df, subs=subs, decompress = T)
  rm(subs)
  
  #add template
  df$Template <- "3'U4"
  for(i in 1:length(C)){
    df$Template[df$Subunit==C[i]] <- "3'C4"
  }
  rm(C)
  
  #factor for glm
  df$Template <- factor(df$Template, levels=c("3'U4", "3'C4"))
  df$Subunit <- factor(df$Subunit, levels=c(subNames))
  rm(subNames)
  
  df$strain <- strain
  
  #make it look pretty
  df <- df[,c(
    grep('series', colnames(df)), 
    grep('rounds', colnames(df)),
    grep('strain', colnames(df)),
    grep('Subunit', colnames(df)),
    grep('Template', colnames(df)),
    grep('Seq_Trim_R', colnames(df))
    )
    ]
  return(df)
}

#setwd("") #location of the data if not default directory
names <- c('pr8', 'hk', 'wsn', 'bri')
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
cluster <- length(names)
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar% 
  dfPreProcess(read.csv(paste(toupper(names[j]), '_All_Match.csv', sep=''), stringsAsFactors = F), strains[j], NA, NA)
stopCluster(cl)

setwd("./Overview")

#vectors for the total number of mRNA, the number in no realignment (no) and realigned (par)
reads <- vector(mode='numeric', length=length(strains))
no <- vector(mode='numeric', length=length(strains))
par <- vector(mode='numeric', length=length(strains))

for(i in 1:length(strains)){
  reads[i] <- nrow(dfs[[i]])
  no[i] <- nrow(dfs[[i]][dfs[[i]]$series=='No_Realignment',])
  par[i] <- nrow(dfs[[i]][dfs[[i]]$series=='Realigned',])
}

#par frequency data.frame
parRates <- matrix(ncol=length(strains), nrow=2)
parRates[1,] <- no/reads * 100
parRates[2,] <- par/reads * 100

#names
row.names(parRates) <- c('No Realignment', 'Realigned')
colnames(parRates) <- strains

write.csv(parRates, 'PAR rates.csv')

##### Stuff with Rounds of PAR
#extract the highest number of rounds in the 4 data sets
maxRounds <- max(dfs[[1]]$rounds)
for(i in 2:length(strains)){
  maxRounds <- max(maxRounds, dfs[[i]]$rounds)
}

#extract the number of mRNA that undergo the jth number of rounds
parMulti <- matrix(ncol=length(strains), nrow=maxRounds)
for(i in 1:ncol(parMulti)){
  for(j in 1:nrow(parMulti)){
    parMulti[j,i] <- nrow(dfs[[i]][dfs[[i]]$rounds==j,])/par[i] * 100
  }
}
row.names(parMulti) <- paste('Round', 1:nrow(parMulti))
colnames(parMulti) <- strains
write.csv(parMulti, 'PAR rounds.csv') #save data

#PAR survival rates (total mRNA that underfo the ith round of more)
par <- data.frame(parRates)
for(i in 2:nrow(parMulti)){
  if(i == nrow(parMulti)){
    par[nrow(par)+1,] <- as.numeric(parMulti[i:nrow(parMulti),])*parRates[2,]/100
  } else {
    par[nrow(par)+1,] <- as.numeric(colSums(parMulti[i:nrow(parMulti),]))*parRates[2,]/100
  }
}

write.csv(par, 'PAR rounds survival.csv') #save

parMulti <- par[2:nrow(par),] #drops the data point of mRNA that undergo 0 rounds of PAR

#Log transformed data
roundsPlotdf <- matrix(ncol=3, nrow=length(strains)*(maxRounds))
for(i in 1: length(strains)){
  roundsPlotdf[((i-1)*(maxRounds)+1):((i)*(maxRounds)),1] <- strains[i]
  roundsPlotdf[((i-1)*(maxRounds)+1):((i)*(maxRounds)),2] <- 1:maxRounds
  roundsPlotdf[((i-1)*(maxRounds)+1):((i)*(maxRounds)),3] <- parMulti[,i]
}
roundsPlotdf <- data.frame(strain=roundsPlotdf[,1], round=as.numeric(roundsPlotdf[,2]), percent=log(as.numeric(roundsPlotdf[,3])))

#drop NA and -Inf
roundsPlotdf[,3][is.infinite(roundsPlotdf[,3])] <- NA
roundsPlotdf <- roundsPlotdf[!is.na(roundsPlotdf[,3]),]

#factor for plot
roundsPlotdf[,1] <- factor(roundsPlotdf[,1], levels = strains)

#get R^2
r2val <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  r2val[i] <- summary(lm(data = roundsPlotdf[roundsPlotdf$strain==strains[i],], round ~ percent))$r.squared
}

#generate rounds correlation plot on log scale
roundsPlot <- ggplot(roundsPlotdf, aes(x=round, y=percent, group=strain, colour=strain)) +
  geom_point(size=0.5) +
  geom_smooth(method='lm', se=F, formula = y~x, size=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(1, max(roundsPlotdf$round), by=1)) +
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
  labs(x='Rounds of PAR', y='Realigned Primers (ln(%))')

#steal legend
roundsLegend <- get_legend(roundsPlot)

#remove legend
roundsPlot <- roundsPlot + theme(legend.position = 'none')

#save
save_plot('Rounds_Log_Survival.png', plot=plot_grid(roundsPlot, roundsLegend, rel_widths = c(2.5,1)), base_height = 2, base_width = 3.5, dpi = 600)

#get PAR rate after round 1 from slope of lines
lmPARrates <- vector(mode = 'numeric', length=length(strains))
for(i in 1:length(strains)){
  lmPARrates[i] <- exp(predict(lm(data = roundsPlotdf[roundsPlotdf$strain==strains[i],], percent ~ round), newdata = data.frame(round=2))) / exp(predict(lm(data = roundsPlotdf[roundsPlotdf$strain==strains[i],], percent ~ round), newdata = data.frame(round=1))) *100
}
