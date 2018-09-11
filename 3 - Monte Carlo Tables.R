#setwd() #where the data is if not the default directory

#get PAR frequencies
setwd('./Overview') #requires 1 - PAR Frequencies and Rounds
parav <- read.csv('PAR rates.csv', stringsAsFactors = F, row.names = 1)[2,]
pars <- as.numeric(parav[1,])

#get data
setwd('..')
setwd('./G+1/p-values')
pr <- read.csv('Puerto Rico global obtained z and p.csv', row.names = 1, stringsAsFactors = F)
hk <- read.csv('Hong Kong global obtained z and p.csv', row.names = 1, stringsAsFactors = F)
ws <- read.csv('WSN global obtained z and p.csv', row.names = 1, stringsAsFactors = F)
br <- read.csv('Brisbane global obtained z and p.csv', row.names = 1, stringsAsFactors = F)

generateTable <- function(pr, hk, ws, br, position, nucleotide, pars){
  #extract position from Dangling, End, Upstream (-2, -1, +1)
  pr <- pr[,grep(position, colnames(pr))]
  hk <- hk[,grep(position, colnames(hk))]
  ws <- ws[,grep(position, colnames(ws))]
  br <- br[,grep(position, colnames(br))]
  
  #get nucleotide
  nucleotide <- paste('_', nucleotide, sep='') #ensure correct function of grep
  pr <- pr[,grep(nucleotide, colnames(pr))]
  hk <- hk[,grep(nucleotide, colnames(hk))]
  ws <- ws[,grep(nucleotide, colnames(ws))]
  br <- br[,grep(nucleotide, colnames(br))]
  
  fixData <- function(df, nucleotide){
    df <- df[c(1,2,3,9),] #delimits to composition, simulated mean, simulated sd, p-score
    df[2,] <- df[2,1] #simulated mean on all
    df[3,] <- df[3,1] #simulated sd on all
    colnames(df) <- c('All', 'No', 'PAR')
    row.names(df) <- c(paste(gsub('_', '', nucleotide),'%', sep=''), 'Sim%', 'Sim sd', 'p-value')
    return(df)
  }
  #obtain data
  pr <- fixData(pr, nucleotide)
  hk <- fixData(hk, nucleotide)
  ws <- fixData(ws, nucleotide)
  br <- fixData(br, nucleotide)
  
  #get PAR rate
  pr[nrow(pr)+1,] <- (pr[1,3]/100)*(pars[1]/100)/(pr[1,1]/100)*100
  hk[nrow(hk)+1,] <- (hk[1,3]/100)*(pars[2]/100)/(hk[1,1]/100)*100
  ws[nrow(ws)+1,] <- (ws[1,3]/100)*(pars[3]/100)/(ws[1,1]/100)*100
  br[nrow(br)+1,] <- (br[1,3]/100)*(pars[4]/100)/(br[1,1]/100)*100
  
  #merge data
  df <- cbind(pr, hk, ws, br)
  row.names(df)[5] <- 'PAR%'
  return(df)
}

position <- 'Upstream'
nucleotide <- 'G'
GUpstream <- generateTable(pr, hk, ws, br, position, nucleotide, pars)
write.csv(GUpstream, 'G+1.csv')

position <- 'Dangling'
nucleotide <- 'C'
CDangling <- generateTable(pr, hk, ws, br, position, nucleotide, pars)
write.csv(CDangling, 'C-2.csv')

position <- 'End'
nucleotide <- 'A'
AEnd <- generateTable(pr, hk, ws, br, position, nucleotide, pars)
write.csv(AEnd, 'A-1.csv')

position <- 'End'
nucleotide <- 'U'
UEnd <- generateTable(pr, hk, ws, br, position, nucleotide, pars)
write.csv(UEnd, 'U-1.csv')

position <- 'End'
nucleotide <- 'C'
CEnd <- generateTable(pr, hk, ws, br, position, nucleotide, pars)

#AG and CAG preference
#get cleavage sites
setwd('..')
setwd('./Cleavage Sites') #requires 3 - All Cleavage Sites
pzp <- read.csv('Puerto Rico global obtained z and p.csv', row.names = 1, stringsAsFactors = F)
hzp <- read.csv('Hong Kong global obtained z and p.csv', row.names = 1, stringsAsFactors = F)
wzp <- read.csv('WSN global obtained z and p.csv', row.names = 1, stringsAsFactors = F)
bzp <- read.csv('Brisbane global obtained z and p.csv', row.names = 1, stringsAsFactors = F)


csTable <- function(pzp, hzp, wzp, bzp, regex, pars){
  x <- as.data.frame(
    cbind(pzp[,grep(regex, colnames(pzp))[1]], 
          hzp[,grep(regex, colnames(hzp))[1]],
          wzp[,grep(regex, colnames(wzp))[1]],
          bzp[,grep(regex, colnames(bzp))[1]]))
  colnames(x) <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
  row.names(x) <- row.names(pzp)
  
  x[nrow(x)+1,] <- (x[3,]/100)*(pars/100)/(x[1,]/100)*100
  row.names(x)[nrow(x)] <- 'PAR%'
  
  nts <- vector(mode='numeric', length=0)
  sim <- vector(mode='numeric', length=0)
  ssd <- vector(mode='numeric', length=0)
  pvs <- vector(mode='numeric', length=0)
  par <- vector(mode='numeric', length=0)
  
  for(i in 1:4){
    nts <- c(nts, as.numeric(t(x[1:3,1:4])[i,]))
    sim <- c(sim, rep(x[4,i], 3))
    ssd <- c(ssd, rep(x[5,i], 3))
    pvs <- c(pvs, as.numeric(t(x[12:14,1:4])[i,]))
    par <- c(par, rep(x[15,i], 3))
  }
  
  df <- matrix(nrow=5, ncol=length(nts))
  df[1,] <- nts
  df[2,] <- sim
  df[3,] <- ssd
  df[4,] <- pvs
  df[5,] <- par
  
  df <- as.data.frame(df)
  row.names(df) <- c('Seq%', 'Sim%', 'Sim sd', 'p-value', 'PAR%')
  colnames(df) <- rep(c('All', 'No', 'PAR'),4)
  
  return(df)
}

regex <- '[^C]AG.'
AG <- csTable(pzp, hzp, wzp, bzp, regex, pars)
write.csv(AG, 'AG dinucleotide.csv')

regex <- 'CAG.'
CAG <- csTable(pzp, hzp, wzp, bzp, regex, pars)
write.csv(CAG, 'CAG trinucleotide.csv')