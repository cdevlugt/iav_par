#one of the first scripts I wrote for my MSc
require(ggplot2)
require(reshape2)
require(cowplot)

keepR <- function(table, col=1, val=0, val2=NA, t="==", chr=FALSE){
  # Keeps Rows (keepR) within a specified range
  #
  # Args:
  #	table:	a data.frame to perform operations on, this data.frame will have only rows
  #			except that meet the specified parameters
  #	col: 	column number which the range of values kept will be applied to
  #	val: 	value to keep, or rnage limiter number 1
  #	val2:	range limiter number 2; will be set as val by default if needed
  #	t:  type of operation to be performed to keep rows where this value is TRUE
  # chr: spefies if the values kept are to be relted to string length
  #
  # Returns:
  #	table:	the table passed following removal of rows where "t" argued opperation is false
  if(!chr) {
    if(t=="=="){
      table <- table[table[,col]==val,]
    } else if(t=="!=") {
      table <- table[table[,col]!=val,]
    } else if(t=="<") {
      table <- table[table[,col]<val,]
    } else if(t=="<=") {
      table <- table[table[,col]<=val,]
    } else if(t==">") {
      table <- table[table[,col]>val,]
    } else if(t==">=") {
      table <- table[table[,col]>=val,]
    } else if(t=="na") {
      table <- table[is.na(table[,col]),]
    } else if(t=="!na") {
      table <- table[!(is.na(table[,col])),]
    } else if(t=="<>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]<val) & (table[,col]>val2),]
    } else if(t=="><") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]>val) | (table[,col]<val2),]
    } else if(t=="<=>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]<=val) & (table[,col]>=val2),]
    } else if(t==">=<") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]>=val) | (table[,col]<=val2),]
    } else {
      print("invalid type (t)")
    }
  } else if(chr) {
    if(t=="=="){
      table <- table[nchar(table[,col])==val,]
    } else if(t=="!=") {
      table <- table[nchar(table[,col])!=val,]
    } else if(t=="<") {
      table <- table[nchar(table[,col])<val,]
    } else if(t=="<=") {
      table <- table[nchar(table[,col])<=val,]
    } else if(t==">") {
      table <- table[nchar(table[,col])>val,]
    } else if(t==">=") {
      table <- table[nchar(table[,col])>=val,]
    } else if(t=="<>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])<val) & (nchar(table[,col])>val2),]
    } else if(t=="><") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])>val) | (nchar(table[,col])<val2),]
    } else if(t=="<=>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])<=val) & (nchar(table[,col])>=val2),]
    } else if(t==">=<") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])>=val) | (nchar(table[,col])<=val2),]
    } else {
      print("invalid type (t)")
    }
  }
  return(table)
}
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
snlc <- function(df, sub=1, Length=25, Trim_Len=26, Len_Adjust=27, min=5, max=25, freq=TRUE){
  # Subinit nt_length comparison 
  #
  # Args:
  #	df:	a data.frame to perform operations on
  # sub:  subunit data to return
  # Length: column containing length of full sequence
  # Trim_Len: column containing the length of the fully trimmed sequence
  # Len_Adjust: column containing length adjusted
  #	min: 	length minimum value to tabulate from
  #	max:	length maximum value to tablulate data to
  #	freq:	flag to determine if the data is returned as freqency or as a count
  #
  # Returns:
  #	nt_df:	a df contining the length of snatched and matched sequences within
  #				the range specified by min and max for a given IAV subunit
  
  df <- df[,c(sub, Length, Trim_Len, Len_Adjust)]
  df <- subunitDecompress(df, sub=1, drop=F)
  
  dfStuff <<- df
  
  tlen <- data.frame(table(df[,2]))	#captured sequence (all)
  tlen[,1] <- as.numeric(as.character(tlen[,1]))	#adjusts names back to numeric type
  
  mlen <- data.frame(table(df[,3]))	#matched sequence (all)
  mlen[,1] <- as.numeric(as.character(mlen[,1]))
  
  nal <- data.frame(table(keepR(df, col=4, val=0, t="==")[,2]))	#sequences which were not trimmed to get a match
  nal[,1] <- as.numeric(as.character(nal[,1]))
  
  df <- keepR(df, col=4, val=0, t=">")	#removes rows where there is no realignment
  
  prl <- data.frame(table(df[,2]))	#final length of (after addition) sequences which were matched after atleast 1 round of removal
  prl[,1] <- as.numeric(as.character(prl[,1]))
  
  snl <- data.frame(table(df[,3]))	#length of first match after rounds of removal
  snl[,1] <- as.numeric(as.character(snl[,1]))
  
  add <- data.frame(table(df[,4]))	#number of nt added
  add[,1] <- as.numeric(as.character(add[,1]))
  
  #combines all of those dfs into an intermediate df (int_df)
  #will throw warnings, but it does its job correctly
  int_df <- merge(tlen, mlen, by="Var1", all="TRUE")
  int_df <- merge(int_df, nal, by="Var1", all="TRUE")
  int_df <- merge(int_df, prl, by="Var1", all="TRUE")
  int_df <- merge(int_df, snl, by="Var1", all="TRUE")
  int_df <- merge(int_df, add, by="Var1", all="TRUE")
  int_df[,2:7][is.na(int_df[,2:7])] <- 0
  
  if(freq){	#converts values to frequency if that is specified
    int_df[,2:7] <- int_df[,2:7]/ sum(int_df[,2]) * 100
  }
  
  #creates the df to be returned
  nt_df <- data.frame(Length=seq(min, max), Captured_Sequence=0, Maximally_Trimmed_Sequence=0, No_Realignment=0, Realigned_Sequence=0, Trimmed_Sequence=0, nt_Added=0, stringsAsFactors=FALSE)
  colnames(nt_df) <- c("Length", "Captured", "Maximally Trimmed", "No_Realignment", "Realigned", "Trimmed", "nt added")
  
  for(i in min:max){	#adds data to the df between the specified min and max
    for(j in 1:nrow(int_df)){
      if(i == int_df[j,1]){
        nt_df[(i-min+1),2:7] <- nt_df[(i-min+1),2:7] + int_df[j,2:7]
      }
    }
    nt_df[is.na(nt_df)] <- 0	#if a row is not found, this will get rid of NA errors
  }
  return(nt_df)
}
lengthPlot <- function(df, min=8, max=17, coln="Length"){ 
  # ggplot in a function with ggsave() built into it. accepts a table that has been created using fourColMelt 
  #
  # Args:
  #	df:	a data.frame generated by sMelt; will not work if another table is passed
  #	min:	minimum x-value (minimum nt_length)
  #	max: 	maximum x-value (maximum nt_length)
  # coln: column to melt data to
  #
  #return:
  # graph:  a graph of the passed data
  
  df <- df[,c(1,4,5,6)]	#trims table to columns to analyze
  df[,3:4] <- df[,3:4] * -1	#makes data go to the left of the axis
  df <- melt(df, id=coln)	#melts data
  df$variable <- factor(df$variable, levels = c('No_Realignment', 'Trimmed', 'Realigned'))
  
  graph <- ggplot(data=df, aes(group=variable, x=Length, y=value, fill=variable)) + 
    geom_density(alpha=.9, stat='identity', position='identity') +
    ylab("Frequency (%)") + 
    scale_y_continuous(breaks=seq(-20, 60, by=10), limits=c(-20, 60), labels=c("20", "10", "0", "10", "20", "30", "40", "50", "60")) +
    scale_x_continuous(breaks=seq(min,max, by=1), limits=c(min,max)) +
    xlab("Length (nt)") +
    coord_flip() +
    labs(fill="Sequence Type") +
    theme_bw() +
    scale_fill_manual(name='legend', breaks=c('No_Realignment', 'Trimmed', 'Realigned'), 
                      values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Trimmed=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65), Realigned=hcl(h=seq(15,375, length=(3+1))[3], c=100, l=65)),
                      labels=c(No_Realignment='No Realignment', Trimmed='Pre-Realignment', Realigned='Post-Realignment'))
  
  return(graph)
}
graphLengths <- function(df, subs=c(1:8), min=8, max=17, coln="Length", global=T, subName='Subunit'){
  #Analyzes data in the program and generates boxplots
  #
  #args:
  # df: a data frame containing data
  # subs: columns containing IAV subuunits
  # min:  minimum y-unit
  # max:  maximum x-unit
  # coln: column name to melt around
  # global: should the global values also be plotted
  #
  #return:
  # graph:  a cowplot of density graphs for all selected subunits with titles
  
  #add the g-comp
  df$Trim_Len[startsWith(df$NT_coupure, "G")] <- df$Trim_Len[startsWith(df$NT_coupure, "G")] +1
  df$Length[startsWith(df$NT_coupure, "G") & df$rounds==0] <- df$Length[startsWith(df$NT_coupure, "G") & df$rounds==0] +1
  df$Len_Adjust <- df$Length - df$Trim_Len
  
  if(global){
    df$All <- rowSums(df[,c(subs)])
    subs <- c(ncol(df), subs)
  }
  
  letterNames <- LETTERS[1:length(subs)]
  graphNames <- vector(mode="character", length=0)
  
  #generate graphs for all subunits and global
  for(i in 1:length(subs)){
    print(subs[i])
    if(global && i==1){
      graphNames[i] <- paste(colnames(df)[subs[i]],"combined", sep="_")
      title <- paste("All ", subName,"s",sep="")
    } else {
      if(grepl("\\.", colnames(df)[subs[i]], perl=T)){
        graphNames[i] <- paste(substr(colnames(df)[subs[i]],1,regexpr("\\.",colnames(df)[subs[i]], perl=T)-1),"combined", sep="_")
        title <- paste(substr(colnames(df)[subs[i]],1,regexpr("\\.",colnames(df)[subs[i]], perl=T)-1),subName, sep=" ")
      } else {
        graphNames[i] <- paste(colnames(df)[subs[i]],"combined", sep="_")
        title <- paste(colnames(df)[subs[i]],subName, sep=" ")
      }
    }
    title <- ggdraw() + draw_label(title, fontface='bold')
    if(i==1){
      assign(graphNames[i], lengthPlot(snlc(df, sub=subs[i], Length = grep("Length", colnames(df)), Trim_Len=grep("Trim_Len", colnames(df)), Len_Adjust=grep("Len_Adjust", colnames(df)), min=min, max=max, freq=T), min=min, max=max, coln=coln))
      plotLegend <- get_legend(get(graphNames[i]) + theme(legend.position='left'))
      assign(graphNames[i], plot_grid(title, get(graphNames[i]) + theme(legend.position='none'), labels=c("",letterNames[i]), ncol=1, rel_heights=c(0.1, 1)))
    } else {
      assign(graphNames[i], plot_grid(title, lengthPlot(snlc(df, sub=subs[i], Length = grep("Length", colnames(df)), Trim_Len=grep("Trim_Len", colnames(df)), Len_Adjust=grep("Len_Adjust", colnames(df)), min=min, max=max, freq=T), min=min, max=max, coln=coln) + theme(legend.position='none'), labels=c("",letterNames[i]), ncol=1, rel_heights=c(0.1, 1)))
    }
  }
  
  graph <- plot_grid(get(graphNames[1]), get(graphNames[2]),get(graphNames[3]), get(graphNames[4]), plotLegend, ncol=5, nrow=1, rel_widths = c(1,1,1,1,.3))
  
  if(length(graphNames)>4){
    for(i in 2:ceiling(length(graphNames)/4)){
      graph <- plot_grid(graph, plot_grid(get(graphNames[(4*(i-1)+1)]),get(graphNames[(4*(i-1)+2)]),get(graphNames[(4*(i-1)+3)]), get(graphNames[(4*(i-1)+4)]), plotLegend + theme(legend.position='none'),ncol=5, nrow=1, rel_widths = c(1,1,1,1,.3)), ncol=1, rel_heights = c((i-1),1))
    }
  }
  return(graph)
}

#data directory
setwd("")
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

#you will need to create a Length folder in this directory
setwd("./Length")

#genome segments length
pr8plot <- graphLengths(pr8, global=F)
save_plot("PR8_Length_Subunits_G.png",pr8plot, base_height = 9, base_width = 26*.75, dpi=900)

hkplot <- graphLengths(hk, global=F)
save_plot("hk_Length_Subunits_G.png",hkplot, base_height = 9, base_width = 26*.75, dpi=900)

briplot <- graphLengths(bri, global=F)
save_plot("bri_Length_Subunits_G.png",briplot, base_height = 9, base_width = 26*.75, dpi=900)

wsnplot <- graphLengths(wsn, global=F)
save_plot("wsn_Length_Subunits_G.png",wsnplot, base_height = 9, base_width = 26*.75, dpi=900)

#length on the strain level. columns are created for all of the strains, and then the data frames are joined
pr8[,(ncol(pr8)+1)] <- rowSums(pr8[,c(1:8)])
pr8[,(ncol(pr8)+1)] <- 0
pr8[,(ncol(pr8)+1)] <- 0
pr8[,(ncol(pr8)+1)] <- 0
colnames(pr8)[(ncol(pr8)-3):ncol(pr8)] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")
pr8 <- pr8[,-c(1:8, grep("_Trim_R", colnames(pr8)))]
pr8 <- pr8[,c((ncol(pr8)-3):ncol(pr8), 1:(ncol(pr8)-4))]

hk[,(ncol(hk)+1)] <- 0
hk[,(ncol(hk)+1)] <- rowSums(hk[,c(1:8)])
hk[,(ncol(hk)+1)] <- 0
hk[,(ncol(hk)+1)] <- 0
colnames(hk)[(ncol(hk)-3):ncol(hk)] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")
hk <- hk[,-c(1:8, grep("_Trim_R", colnames(hk)))]
hk <- hk[,c((ncol(hk)-3):ncol(hk), 1:(ncol(hk)-4))]

wsn[,(ncol(wsn)+1)] <- 0
wsn[,(ncol(wsn)+1)] <- 0
wsn[,(ncol(wsn)+1)] <- rowSums(wsn[,c(1:8)])
wsn[,(ncol(wsn)+1)] <- 0
colnames(wsn)[(ncol(wsn)-3):ncol(wsn)] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")
wsn <- wsn[,-c(1:8, grep("_Trim_R", colnames(wsn)))]
wsn <- wsn[,c((ncol(wsn)-3):ncol(wsn), 1:(ncol(wsn)-4))]

bri[,(ncol(bri)+1)] <- 0
bri[,(ncol(bri)+1)] <- 0
bri[,(ncol(bri)+1)] <- 0
bri[,(ncol(bri)+1)] <- rowSums(bri[,c(1:8)])
colnames(bri)[(ncol(bri)-3):ncol(bri)] <- c("Puerto Rico", "Hong Kong", "WSN", "Brisbane")
bri <- bri[,-c(1:8, grep("_Trim_R", colnames(bri)))]
bri <- bri[,c((ncol(bri)-3):ncol(bri), 1:(ncol(bri)-4))]

global <- rbind(pr8, hk, wsn, bri)
rm(pr8, hk, wsn, bri)

globalPlot <- graphLengths(global, subs=c(1:4), global=F, subName = "Strain")
save_plot("Global_Length_G.png",globalPlot, base_height = 4.5, base_width = 26*.75, dpi=900)
