require(ggplot2)
require(cowplot)
require(ggseqlogo)

require(compiler)
enableJIT(3)

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
cleavageSite <- function(df, range=3, subs=c(1:8), seq=11, trim=15, downstream=16, seperate=T, decompress=T, name="Subunit"){
  
  df <- df[,c(subs, seq, trim, downstream)]
  
  subs <- 1:length(subs)
  seq <- length(subs)+1
  trim <- seq + 1
  downstream <- trim+1
  
  df$TrimFlag <- FALSE
  df$TrimFlag[df[,seq]!=df[,trim]] <- TRUE
  
  df$Cleave <- paste(substr(df[,trim], nchar(df[,trim])-(range-1), nchar(df[,trim])), substr(df[,downstream], 1,range),sep="")
  df$Length <- nchar(df[,trim])
  
  df <- df[,c(subs, ncol(df)-1, ncol(df), ncol(df)-2)]
  rm(seq, trim, downstream)
  
  if(seperate){
    df <- subunitSeperate(df, subs=subs, decompress = decompress, name=name)
  }
  
  return(df)
}
splitLogo <- function(x, title=NA){
  
  if(length(x)==0){
    return(NULL)
  }
  
  splitPoint <- nchar(x[1])/2
  
  sl <- ggplot()+
    geom_logo(x, method='prob') +
    theme_logo() +
    annotate('rect', xmin=.5, xmax=splitPoint+.5, ymin=0, ymax=1, col='#000000', size=1.5, fill="transparent") +
    annotate('rect', xmin=splitPoint+.5, xmax=2*splitPoint+.5, ymin=0, ymax=1, col='#000000', size=1.5, fill="transparent") +
    annotate('text', x=(splitPoint+1)/2, y=1.1, label="Sequence End", size=3) +
    annotate('text', x=(splitPoint+1)/2+splitPoint, y=1.1, label="Downstream NT", size=3) +
    scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
    scale_x_continuous(breaks=c(1:splitPoint,(splitPoint+1):(2*splitPoint)), labels = c((-1*splitPoint):-1, 1:splitPoint)) +
    ylab("Frequency") +
    xlab("Position") +
    theme(axis.title = element_text(size=8), axis.text = element_text(size=6))
  
  if(!is.na(title)){
    sl <- sl + ggtitle(title) +
      theme(plot.title = element_text(face='bold', size=10, hjust=0.5))
  }
  
  return(sl)
}
realignmentFrequency <- function(df=NULL, freq=T, name=''){
  
  if(is.null(nrow(df))){
    re <- data.frame(series=c('No_Realignment', 'Trimmed'), x=1, y=1)
    
    reGraph <- ggplot() +
      geom_bar(data=re, aes(group=series, y=y, x=x, fill=series), stat='identity') + 
      theme_bw() +
      scale_fill_manual(name='Series', breaks=c('No_Realignment', 'Trimmed'), 
                        values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Trimmed=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65)),
                        labels=c(No_Realignment='No Realignment', Trimmed='Realigned'))
    return(get_legend(reGraph))
  } else if(nrow(df)==0) {
    return(NULL)
  }
  
  #receives dfTempY from graphCleavageSite
  
  df$end <- substr(df$Cleave, 4, 4)
  df$end <- gsub('[^G]', '', df$end)
  
  df$series <- 'No Realignment'
  df$series[df$TrimFlag==T] <- 'Trimmed'
  
  #generate dataframe with realignment frequencies
  re <- data.frame(
    series = c(rep('No_Realignment',2), rep('Trimmed',2), rep('all', 2)),
    end = c('G', ''),
    y=0, 
    x=c(1:2),
    stringsAsFactors = F
  )
  
  #get values
  for(i in 1:2){
    re[i,3] <- nrow(df[df$series=='No Realignment' & df$end==re[i,2],])
    re[(i+2),3] <- nrow(df[df$series=='Trimmed' & df$end==re[(i+2),2],])
    re[(i+4),3] <- nrow(df[df$end==re[(i+4),2],])
  }
  
  #write fold changes
  output <- data.frame(x1 = re[,2][re[,1]=='No_Realignment'],
                       x2 = re[,2][re[,1]=='Trimmed'],
                       NoR = re$y[re[,1]=='No_Realignment'],
                       Trim = re$y[re[,1]=='Trimmed'],
                       all = re$y[re[,1]=='all'],
                       PAR = re$y[re[,1]=='Trimmed']/re$y[re[,1]=='all']*100
  )
  write.csv(output, paste(name, '.csv', sep=''))
  rm(output)
  gc()
  
  #get frequencies
  if(freq){
    re[1:2,3] <- re[1:2,3] / re[5:6,3] *100
    re[3:4, 3] <- re[3:4,3] / re[5:6,3] *100
  }
  
  if(freq){
    #get global realignment frequency
    reFreq <- nrow(df[df$series=='Trimmed',])/nrow(df) * 100
  } else {
    #global realignment frequency as count
    reFreq <- data.frame(y=(re[17:24,3] * nrow(df[df$series=='Trimmed',])/nrow(df)), x=0)
    reFreq$x <- 1:nrow(reFreq)
  }
  
  re <- re[-c(5:6),] #drop useless data
  
  #factor for order
  re[,1] <- factor(re[,1], levels=c('No_Realignment', 'Trimmed'))
  
  reGraph <- ggplot() +
    geom_bar(data=re, aes(group=series, y=y, x=x, fill=series), stat='identity') + 
    theme_bw() +
    scale_x_continuous(breaks=c(1:2), labels = c('G', 'Other')) +
    xlab("") +
    scale_fill_manual(breaks=c('No_Realignment', 'Trimmed'), 
                      values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Trimmed=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65)),
                      labels=c(No_Realignment='No Realignment', Trimmed='Realigned'), guide=F) 
  
  if(freq){
    reGraph <- reGraph +
      geom_hline(yintercept = reFreq) +
      ylab('Frequency (%)')
  } else {
    reGraph <- reGraph +
      stat_summary(fun.y='mean', data=reFreq, colour="#000000", geom="errorbar", aes(x=x, y=y,ymax=..y.., ymin=..y..), width=1, linetype="solid", show.legend=F) +
      ylab('Count')
  }
  
  return(reGraph)
}

#setwd("") #data directory
  hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
  pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
  wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
  bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

graphCleavageSite <- function(df, by='subunit', lengths=NA, range=3, subs=c(1:8), seq=11, trim=15, downstream=16, name="Subunit", strain=NA){
  
  
  #sets up analysis loops 
  if(by!='subunit' && by!='length' && by!='both'){
    warning("invalid 'by' selection; options are 'subunit', 'length', or 'both'; set to 'subunit' by default")
  } 
  if(by=='both'){
    x <- length(subs)
    if(is.na(lengths)){
      lengths <- min(nchar(df[,trim])):max(nchar(df[,trim]))
    }
    y <- length(lengths)
    letterNames <- vector(mode ='character', length=4*y)
  } else if(by=='length') {
    x <- 1
    if(is.na(lengths)){
      lengths <- min(nchar(df[,trim])):max(nchar(df[,trim]))
    }
    y <- length(lengths)
    letterNames <- vector(mode ='character', length=4*y)
  } else { #else if 'subunit'
    x <- length(subs)
    y <- 1
    letterNames <- vector(mode ='character', length=4*x)
  }
  
  for(i in 1:ceiling(length(letterNames)/26)){ #constructs letterNames
    if(i==1){
      firstLetter <- ""
    } else {
      firstLetter <- LETTERS[i-1]
    }
    for(j in ((i-1)*26+1):(i*26)){
      if(j>length(letterNames)){
        break
      } else {
        letterNames[j] <- paste(firstLetter, LETTERS[(j-26*(i-1))], sep="")
      }
    }
  }
  rm(firstLetter)
  
  subNames <- colnames(df)[subs]
  df <- cleavageSite(df, range=range, subs=subs, seq=seq, trim=trim, downstream = downstream, seperate=T, decompress=T, name = name) #processes df
  #Subunit|Cleave|Length|TrimFlag
  #Strain|Cleave|Length|TrimFlag
  
  #generate output
  if(by!='both'){
    if(by=='subunit' && name!='Subunit'){
      by = 'Strain'
    }
    outputVar <- unique(df[,grep(by, colnames(df), ignore.case = T)])
    if(is.integer(outputVar)){
      outputVar <- outputVar[order(outputVar)]
    }
    output <- data.frame(
      target = outputVar,
      No_Realignment = 0,
      Realigned = 0,
      stringsAsFactors = F
    )
    rm(outputVar)
    gc()
    
    for(i in 1:nrow(output)){ #this gets all rows
      output[i,2] <- nrow(df[df[,grep(by, colnames(df), ignore.case = T)]==output$target[i] & df$TrimFlag==F & grepl('...G..', df$Cleave),])/nrow(df[df[,grep(by, colnames(df), ignore.case = T)]==output$target[i] & df$TrimFlag==F,]) *100
      output[i,3] <- nrow(df[df[,grep(by, colnames(df), ignore.case = T)]==output$target[i] & df$TrimFlag==T & grepl('...G..', df$Cleave),])/nrow(df[df[,grep(by, colnames(df), ignore.case = T)]==output$target[i] & df$TrimFlag==T,]) *100
    }
    write.csv(output, paste(strain, " by ", by, '.csv', sep=''), row.names = F)
    rm(output)
    
  } else { #case of both
    subV <- unique(df$Subunit)
    lenV <- unique(df$Length)
    lenV <- lenV[order(lenV)]
    
    for(j in 1:length(subV)){
      output <- data.frame(
        target = lenV,
        No_Realignment = 0,
        Realigned = 0,
        stringsAsFactors = F
      )
      
      for(i in 1:nrow(output)){ #this gets all rows for that subunit
        output[i,2] <- nrow(df[df$Length==output$target[i] & df$TrimFlag==F & grepl('...G..', df$Cleave) & df$Subunit==subV[j],])/nrow(df[df$Length==output$target[i] & df$TrimFlag==F & df$Subunit == subV[j],]) *100
        output[i,3] <- nrow(df[df$Length==output$target[i] & df$TrimFlag==T & grepl('...G..', df$Cleave) & df$Subunit==subV[j],])/nrow(df[df$Length==output$target[i] & df$TrimFlag==T & df$Subunit == subV[j],]) *100
      }
      print(paste(strain, " ", gsub('\\.', '', subV[j]), " subunit by length", '.csv', sep=''))
      write.csv(output, paste(strain, " ", gsub('\\.', '', subV[j]), " subunit by length", '.csv', sep=''), row.names=F)
      rm(output)
      gc()
    }
    rm(subV, lenV)
  }
  gc()
  
  if(is.na(strain)){
    strain = name
  }
  
  titles <- vector(mode='character', length=(x*y))
  if(x>1 && y>1){
    #both
    for(i in 1:x){
      for(j in 1:y){
        titles[(i-1)*y+j] <- paste(subNames[i]," ", name, " - ", lengths[j], " Nucleotides",sep="")
      }
    }
    fName <- paste(strain,"by","Length_by_Subunit",subNames, name, sep="_")
  } else if(x>1) {
    #subunit
    titles <- paste(subNames, name, sep=" ")
    fName <- paste(strain,"by",by, sep="_")
  } else if(y>1) {
    #length
    titles <- paste(strain, " - ", lengths, " Nucleotides", sep="")
    fName <- paste(strain,"by",by, sep="_")
  } else {
    titles <- strain
    fName <- paste(strain,"by",by, sep="_")
  }
  titles <- gsub("\\.", "", titles)
  fName <- gsub("\\.", "", fName)
  fName <- paste(fName,".png", sep="")
  
  graph <- NA
  saved <- 1
  
  for(i in 1:x){ #subunit handling
    if(x>1){
      dfTempX <- df[grep(subNames[i], df[,1]),]
      count=nrow(dfTempX)
    } else if(x==1 && y==1) {
      dfTempX <- df
      rm(df)
      gc()
      subs <- name
    }
    
    for(j in 1:y){ #length handling
      if(y>1 && x>1){
        dfTempY <- dfTempX[dfTempX[,3]==lengths[j],]
      } else if(y>1){
        dfTempY <- df[df[,3]==lengths[j],]
        count=nrow(df)
      } else {
        dfTempY <- dfTempX
        rm(dfTempX)
        gc()
      }
      title <- ggdraw() + draw_label(titles[(i-1)*y+j], fontface='bold', size = 13)
      print(titles[(i-1)*y+j])
      if(j==1 && is.na(graph)){
        graph <- plot_grid(
          splitLogo(dfTempY$Cleave, title="All"), splitLogo(dfTempY$Cleave[dfTempY$TrimFlag==F], title="No Realignment"), splitLogo(dfTempY$Cleave[dfTempY$TrimFlag==T], title="Realigned"), realignmentFrequency(dfTempY, freq=T, titles[(i-1)*y+j]), labels=c(),nrow=1, ncol=4, rel_widths=c(1,1,1,1)) +
          theme(aspect.ratio = 4.5/22)
        graph <- plot_grid(title, graph, ncol=1, nrow=2, rel_heights = c(.75, 4.5))
        graphNum <- 1
        
      } else {
        graphTemp <- plot_grid(splitLogo(dfTempY$Cleave, title="All"), splitLogo(dfTempY$Cleave[dfTempY$TrimFlag==F], title="No Realignment"), splitLogo(dfTempY$Cleave[dfTempY$TrimFlag==T], title="Realigned"), realignmentFrequency(dfTempY, freq=T, titles[(i-1)*y+j]), nrow=1, ncol=4, rel_widths=c(1,1,1,1)) + 
          theme(aspect.ratio = 4.5/22)
        graphTemp <- plot_grid(title, graphTemp, ncol=1, nrow=2, rel_heights = c(.75, 4.5))
        graph <- plot_grid(graph, graphTemp, ncol=1, nrow=2, rel_heights = c(1*graphNum, 1))
        rm(graphTemp)
        gc()
        graphNum <- graphNum + 1
      }
      
      rm(dfTempY)
      gc()
    }
    
    if(y!=1){ #removes dfTempX if it has not already been removed
      rm(dfTempX)
      gc()
    }
    if(y>1 && x>1){
      save_plot(fName[saved],graph, base_width=11, base_height=2.625*graphNum, dpi=900)
      graph <- NA
      saved <- saved + 1
    } else if(i==x) {
      save_plot(fName[saved],graph, base_width=11, base_height=2.625*graphNum, dpi=900)
      graph <- NA
      saved <- saved + 1
    }
  }
}

#you will need to create a G+1 folder and a folder for each of the 4 strains
setwd("./G+1/Brisbane")
strain="Brisbane"
graphCleavageSite(bri, strain=strain)
graphCleavageSite(bri, by='length', lengths=9:16, strain=strain)
graphCleavageSite(bri, by='both', lengths=9:16, strain=strain)

setwd('..')
setwd("./Hong Kong")
strain="Hong Kong"
graphCleavageSite(hk, strain=strain)
graphCleavageSite(hk, by='length', lengths=9:16, strain=strain)
graphCleavageSite(hk, by='both', lengths=9:16, strain=strain)

setwd('..')
setwd("./Puerto Rico")
strain="Puerto Rico"
graphCleavageSite(pr8, strain=strain)
graphCleavageSite(pr8, by='length', lengths=9:16, strain=strain)
graphCleavageSite(pr8, by='both', lengths=9:16, strain=strain)

setwd('..')
setwd("./WSN")
strain="WSN"
graphCleavageSite(wsn, strain=strain)
graphCleavageSite(wsn, by='length', lengths=9:16, strain=strain)
graphCleavageSite(wsn, by='both', lengths=9:16, strain=strain)

setwd("..")
#merge strain data and graph all strains
strain="All Strain"
pr8[,ncol(pr8)+1] <- rowSums(pr8[,c(1:8)])
pr8[,ncol(pr8)+1] <- 0
pr8[,ncol(pr8)+1] <- 0
pr8[,ncol(pr8)+1] <- 0
pr8 <- pr8[,c((ncol(pr8)-3):ncol(pr8), 11,15,16)]
colnames(pr8)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

hk[,ncol(hk)+1] <- 0
hk[,ncol(hk)+1] <- rowSums(hk[,c(1:8)])
hk[,ncol(hk)+1] <- 0
hk[,ncol(hk)+1] <- 0
hk <- hk[,c((ncol(hk)-3):ncol(hk), 11,15,16)]
colnames(hk)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

wsn[,ncol(wsn)+1] <- 0
wsn[,ncol(wsn)+1] <- 0
wsn[,ncol(wsn)+1] <- rowSums(wsn[,c(1:8)])
wsn[,ncol(wsn)+1] <- 0
wsn <- wsn[,c((ncol(wsn)-3):ncol(wsn), 11,15,16)]
colnames(wsn)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

bri[,ncol(bri)+1] <- 0
bri[,ncol(bri)+1] <- 0
bri[,ncol(bri)+1] <- 0
bri[,ncol(bri)+1] <- rowSums(bri[,c(1:8)])
bri <- bri[,c((ncol(bri)-3):ncol(bri), 11,15,16)]
colnames(bri)[1:4] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")

All_strain <- rbind(pr8, hk, wsn, bri)
rm(pr8, hk, wsn, bri)
gc()
graphCleavageSite(All_strain,subs=c(1:4), seq=5,trim=6,downstream=7, name="Strain", strain=strain)
