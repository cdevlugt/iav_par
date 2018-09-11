memory.limit(size=(128000+196000)) #phyical RAM + a page file
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
seqComp <- function(df, chars=c('A', 'C', 'G', 'T'), Sequences=10, otherCols=NA, ignore.case=F, Occur=T, Freq=T, All=F, purines=F, strong=F){
  #calculates the relative or absolute portion of a string that is made up of a substring; works to determine the portion of a string made up of
  # RNA or DNA bases/motifs
  #
  #Args:
  # df: dataframe or vector containg the strings
  # chars:  substrings to search for
  # Sequences:  columns containing strings to analyze
  # otherCols:  allows for retention of other data frame columns; you may also cbind results and ignore this feature
  # ignore.case:  ignore case when searching for substrings
  # Occur:  count occurances
  # Freq: count frequency
  # All:  return all data rather than just occurances and/or frequency
  # purines:  determine A and G content of string
  # strong: determine G and C content of string
  #
  #return:
  # cbind(df, dfR): data frame containing sequences and other columns
  # dfR:  dataframe containing sequence composition as specified in Occur and Freq
  
  
  if(!Occur && !Freq){
    stop("Specify either Occur of Freq as true; Occur returns occurances, Freq returns frequency; both may be true")
  }
  
  cols <- Sequences
  Sequences <- as.vector(1:length(Sequences))
  
  if(!is.na(otherCols[1])){
    cols <- c(otherCols, cols)
    otherCols <- 1:length(otherCols)
    Sequences <- as.vector((length(cols)-length(Sequences)+1):length(cols))
  }
  df <- data.frame(df[,cols], stringsAsFactors = F)
  cols <- ncol(df)
  if(Occur){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(chars[j])
        colnames(df)[ncol(df)] <- paste("Abs",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    dfR <- df[,c((cols+1):ncol(df))]
    df <- df[,1:cols]
  }
  if(Freq){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    
    if(Occur){
      dfR <- cbind(dfR, df[,c((cols+1):ncol(df))])
    } else {
      dfR <- df[,c((cols+1):ncol(df))]
    }
    dfR[is.na(dfR)] <-0
    df <- df[,1:cols]
  }
  
  if(All){
    return(cbind(df, dfR))
  } else {
    return(dfR)
  }
}
strMismatch <- function(x=".", mismatch=NA, vect=TRUE, intermediates=TRUE){
  # creates a string containing every possible version of a string with a given
  #   number of mismatches. uses regular expressions for this.
  #
  # Args:
  #   x:  a string or an array of strings
  #   mismatch: number of possible mismatches allowed; forced to run atleast once, so less
  #     than one will still run once; serves to count recursion
  #   vect: specifies if the return should be as a vector rather than a string
  #   intermediates:  logical opperator to return intermediates rather than just the initial string
  #     with n-mutations
  # 
  # Returns:
  #   x:  returns the string if mismatch=0 on first run; doesn't work if a vector is passed
  #	  strMismatch(x.strs, mismatch=mismatch-1, vect=vect, intermediates = intermediates):
  #     Recusion, passing the array of strings to the next itteration, as well as mismatch-1
  #     to count recursion
  #   paste(x.strs, collapse = "|"):
  #     collapses vector into a string
  #   x.strs: an array containing all the generated strings
  
  if(is.na(mismatch)){  #handles passing of a zero mismatch on first run through
    return(x)
  }
  
  x.strs <- vector(mode="character", length=0)
  for(j in 1:length(x)){
    x.chrs <- unlist(strsplit(x[j], split=""))
    x.index <- c(1:nchar(x[j]))
    x.drops <- vector(mode="numeric", length=0)
    for(i in 1:nchar(x[j])){
      if(x.chrs[i]=="."){
        x.drops[length(x.drops)+1] <- i
      } else if(x.chrs[i]=="[") {
        open <- i+1
        close <- regexpr("]", substr(x[j], i, nchar(x[j])))[1]+i-1
        x.drops[(length(x.drops)+1):(length(x.drops)+length(open:close))] <- open:close
      }
    }
    if(length(x.drops)>0){
      x.index <- x.index[-x.drops]
    }
    for(i in 1:length(x.index)){
      if(x.chrs[x.index[i]]=="["){  #handles [stuff] regular expressions
        open <- x.index[i]
        close <- regexpr("]", substr(x[j], x.index[i], nchar(x[j])))[1]+open-1
        x.strs[length(x.strs)+1] <- paste(substr(x[j],1,x.index[i]-1),sub(paste("\\",substr(x[j],open,close),sep=""),".",substr(x[j],x.index[i],nchar(x[j]))), sep="")
      } else {  #handles everything else
        x.str.temp <- x.chrs
        x.str.temp[x.index[i]] <- "."
        x.strs[length(x.strs)+1] <- paste(x.str.temp, collapse="")
      }
    }
  }
  
  #include intermediates and force unique
  if(intermediates){
    x.strs <- unique(c(x, x.strs))
  } else {
    x.strs <- unique(x.strs)
  }
  
  #returns
  if(mismatch<=1){
    if(vect){
      return(x.strs)
    } else {
      return(paste(x.strs, collapse = "|")) #regex
    }
  } else { #recursion
    return(strMismatch(x.strs, mismatch=mismatch-1, vect=vect, intermediates = intermediates))
  }
} #needed to obtain possible formulas

#g included dfpp
dfPreProcess <- function(df, strain='', U=1:5, C=6:8) {
  
  #series handing
  df$series <- 'No_Realignment'
  df$series[df$rounds>=1] <- 'Realigned'
  
  #subunit handling
  subs <- c(U, C)
  subs <- unique(subs)
  subs <- subs[order(subs)]
  subNames <- colnames(df)[subs]
  C <- colnames(df)[C]
  
  #set T as U because mRNA
  df$Trim_Sequence <- gsub('T', 'U', df$Trim_Sequence)
  
  #delimit columns
  df <- df[,c(subs, grep('Trim_Sequence', colnames(df)), grep('NT_coupure', colnames(df)), grep('Trim_Len', colnames(df)), grep('series', colnames(df)))]
  gc()
  
  #G+1
  df$NT_coupure <- substr(df$NT_coupure, 1, 1)
  df$NT_coupure <- gsub('[^G]', "", df$NT_coupure)
  
  #Sequence End
  df$Seq_End <- substr(df$Trim_Sequence, nchar(df$Trim_Sequence), nchar(df$Trim_Sequence))
  df$Seq_End <- paste(df$Seq_End, df$NT_coupure, sep='')
  
  #GCA content
  df$Trim_Sequence[df$NT_coupure=='G'] <- paste(df$Trim_Sequence[df$NT_coupure=='G'], "G", sep='')
  df$GCA_content <- 100 - seqComp(df, chars='U', grep('Trim_Sequence', colnames(df)), Occur=F, Freq=T)
  
  #Length
  df$Trim_Len[df$NT_coupure=='G'] <- df$Trim_Len[df$NT_coupure=='G'] + 1
  
  #organize df
  df <- df[,-(grep('Trim_Sequence', colnames(df)))]
  colnames(df)[(length(subs)+1):(length(subs)+2)] <- c('G_Comp', 'Length')
  
  #since data is G_Comp included, drop G_Comp
  df <- df[,-(grep('G_Comp', colnames(df)))]
  
  #factor for glm
  df$series <- factor(df$series, levels=c('No_Realignment', 'Realigned'))
  
  #length is factored because data is of type int, not double; decimal levels are impractical therefore it must be fitted as a factor
  df$Length <- factor(df$Length, levels=c(9:17))
  df$Seq_End <- factor(df$Seq_End, levels =c('U', 'C', 'G', 'A', 'UG', 'CG', 'GG', 'AG'))
  
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
  #Subunit|Length|series|Seq_End|GCA_content|Template|strain
  df <- df[,c(3,2,4,6,5,1,7)]
  #series|Length|Seq_End|Template|GCA_content|Subunit|strain
  
  return(df)
}
generateFormulas <- function(y, xs){
  #generates all possible formulas including level 2 interactions to predict y
  # from variables xs
  #a formula is redundant if it has only a subset of interaction terms added
  #simply, if a, b, and c are predictors:
  # ~ 1 + a + b + c + a:b + a:c + b:c is best with formulas such as:
  # ~ 1 + a + a:b + a:c + b:c being redundant or objectively inferior
  
  predictorSubset <- function(xs){
    #generates a data.frame of booleans specifying if a variable is on or off
    #can be collapsed into formulas
    #get number of predictor sets
    numx <- length(xs)
    
    #create 'T' array
    tfs <- paste(rep('T', numx), collapse = '')
    tfs <- gsub("\\.", "F", strMismatch(tfs, numx, T, T))
    
    #create matrix
    formFrame <- matrix(nrow=length(tfs), ncol=numx)
    colnames(formFrame) <- xs
    
    #obtain T/F values
    for(i in 1:numx){
      formFrame[,i] <- as.logical(substr(tfs, i, i))
    }
    return(formFrame)
  } #end predictorSubset
  
  lv2Interactions <- function(xs){
    #obtains all permutations of xs
    if(length(xs)==1){
      return(paste(' ~ 1 + ', xs, sep=''))
    } else {
      k = 1
      cxl2 <- vector(mode='character', length=0)
      for(i in 1:(length(xs)-1)){
        for(j in (i+1):(length(xs))){
          cxl2[k] <- paste(xs[i], xs[j], sep=':')
          k = k +1
        }
      }
      
      formFrame <- predictorSubset(cxl2)
      
      #make them a formula
      formulas <- vector(mode='character', length=(nrow(formFrame)-1))
      for(i in 1:length(formulas)){
        formulas[i] <- paste(' ~ 1 + ', paste(xs, collapse = " + "), ' + ', paste(colnames(formFrame)[formFrame[i,]], collapse = " + "), sep='')
      }
      formulas <- c(paste(' ~ 1 + ', paste(xs, collapse = " + "), sep=''), formulas) #formula with no interactions
      
      return(formulas)
    }
  } #end lv2Interactions
  
  #create a table of T/F to pass appropriate variables to lv2Interactions
  formFrame <- predictorSubset(xs)
  
  formulas <- ' ~ 1'
  for(i in 1:(nrow(formFrame)-1)){ #last row is all false
    formulas <- c(formulas, lv2Interactions(colnames(formFrame)[formFrame[i,]])) #passes to lv2Interactions the colnames that have a true value in the ith row
  }
  
  #add y
  formulas <- paste(y, formulas, sep='')
  
  #return the formulas
  return(formulas)
}
slim.glm.aic <- function(name, formula, data, ...) {
  
  cm <- glm(as.formula(paste(deparse(formula))), data=data, ...)
  r2 <- 1 - logLik(cm)/logLik(update(cm, ~ 1)) #McFadden R^2
  rm(data)
  
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  saveRDS(cm, file=paste(name, " ", gsub(':','x',formula),'.rds', sep='')) #save glm
  aic <- cm$aic
  
  rm(cm)
  gc()
  #character vector of formula, aic, and McFadded R2
  return(c(formula, aic, r2))
}

#read bri
#setwd("") #if not default directory
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)
bri <- dfPreProcess(bri, 'Brisbane')
#set bri as df
df <- bri
rm(bri)
df <- df[,c(1:4)]
gc()


setwd("./Effect Modeling")
formulas <- generateFormulas(colnames(df)[1],colnames(df)[-1])
formulas <- c('series ~ 1 + Length + Seq_End + Template + Length:Seq_End + Length:Template + Seq_End:Template + Length:Seq_End:Template', formulas) #manually adding level 3 interaction

#make cluster
cluster <- 30
cl <- makeCluster(cluster)
registerDoParallel(cl)

aics <- foreach(j = 1:length(formulas), .combine = 'rbind', .inorder = F) %dopar%
  slim.glm.aic(name='./Brisbane', formulas[j], df, family=binomial(link='logit'))
stopCluster(cl)

#convert to df
aics <- data.frame(formula = aics[,1], aics=as.numeric(aics[,2]), pseudo_R2=as.numeric(aics[,3]))
write.csv(aics, file='Brisbane aics.csv', row.names = F)

rm(df, aics, formulas)
gc()

#read hk
setwd("..")
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
hk <- dfPreProcess(hk, 'Hong Kong',  U=1:6, C=7:8)
#set hk as df
df <- hk
rm(hk)
df <- df[,c(1:4)]
gc()

setwd("./Effect Modeling")
formulas <- generateFormulas(colnames(df)[1],colnames(df)[-1])
formulas <- c('series ~ 1 + Length + Seq_End + Template + Length:Seq_End + Length:Template + Seq_End:Template + Length:Seq_End:Template', formulas) #manually adding level 3 interaction

#make cluster
cluster <- 2
cl <- makeCluster(cluster)
registerDoParallel(cl)

aics <- foreach(j = 1:length(formulas), .combine = 'rbind', .inorder = F) %dopar%
  slim.glm.aic(name='./Hong Kong', formulas[j], df, family=binomial(link='logit'))
stopCluster(cl)

#convert to df
aics <- data.frame(formula = aics[,1], aics=as.numeric(aics[,2]), pseudo_R2=as.numeric(aics[,3]))
write.csv(aics, file='Hong Kong aics.csv', row.names = F)

rm(df, aics, formulas)
gc()

#read pr8
setwd("..")
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
pr8 <- dfPreProcess(pr8, 'Puerto Rico')
#set pr8 as df
df <- pr8
rm(pr8)
df <- df[,c(1:4)]
gc()

setwd("./Effect Modeling")
formulas <- generateFormulas(colnames(df)[1],colnames(df)[-1])
formulas <- c('series ~ 1 + Length + Seq_End + Template + Length:Seq_End + Length:Template + Seq_End:Template + Length:Seq_End:Template', formulas) #manually adding level 3 interaction

#make cluster
cluster <- 1
cl <- makeCluster(cluster)
registerDoParallel(cl)

aics <- foreach(j = 1:length(formulas), .combine = 'rbind', .inorder = F) %dopar%
  slim.glm.aic(name='./Puerto Rico', formulas[j], df, family=binomial(link='logit'))
stopCluster(cl)

#convert to df
aics <- data.frame(formula = aics[,1], aics=as.numeric(aics[,2]), pseudo_R2=as.numeric(aics[,3]))
write.csv(aics, file='Puerto Rico aics.csv', row.names = F)

#read wsn
setwd("..")
wsn <- read.csv("wsn_All_Match.csv", stringsAsFactors = F)
wsn <- dfPreProcess(wsn, 'WSN')
#set wsn as df
df <- wsn
rm(wsn)
df <- df[,c(1:4)]
gc()

setwd("./Effect Modeling")
formulas <- generateFormulas(colnames(df)[1],colnames(df)[-1])
formulas <- formulas[!grepl(':', formulas)]

#make cluster
cluster <- 1
cl <- makeCluster(cluster)
registerDoParallel(cl)

aics <- foreach(j = 1:length(formulas), .combine = 'rbind', .inorder = F) %dopar%
  slim.glm.aic(name='./WSN', formulas[j], df, family=binomial(link='logit'))
stopCluster(cl)

#convert to df
aics <- data.frame(formula = aics[,1], aics=as.numeric(aics[,2]), pseudo_R2=as.numeric(aics[,3]))
write.csv(aics, file='WSN aics.csv', row.names = F)

rm(df, aics, formulas)
gc()
