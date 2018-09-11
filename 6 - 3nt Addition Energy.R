require(ggplot2)
require(cowplot)
require(compiler)

enableJIT(3)

helixEnergy <- function(x53=NA, y35=NA, offset=0, calc='average'){
  #calculates the free enegry at 37 degrees for a given helix.
  # Turner energy rules (2004) are used for calculations.
  # handles nearest neighbours, terminal mismatches, dangling ends, and single nucleotide mismatches
  #
  #Args:
  # x53 sequence in 5-3
  # y35 sequence in 3-5
  # offset  positive numbers right shift x53, negative numbers right shift y35
  #           allows for dangling ends and starting at a position other than 1
  #
  #Return:
  # df$totalEnergy  vector containing free energies 
  
  #columns are 5'-3' in turner energy rules
  #nearest neghbours Turner 2004 https://rna.urmc.rochester.edu/NNDB/turner04/index.html
  energyTable <- data.frame(AA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.93),
                            AC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.24,NA),
                            AG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.08,NA,-0.55),
                            AU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-1.1,NA,-1.36,NA),
                            CA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.11,NA,NA,NA,NA),
                            CC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-3.26,NA,NA,NA,NA,NA),
                            CG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.36,NA,-1.41,NA,NA,NA,NA),
                            CU=c(NA,NA,NA,NA,NA,NA,NA,NA,-2.08,NA,-2.11,NA,NA,NA,NA,NA),
                            GA=c(NA,NA,NA,NA,NA,NA,NA,-2.35,NA,NA,NA,NA,NA,NA,NA,-1.27),
                            GC=c(NA,NA,NA,NA,NA,NA,-3.42,NA,NA,NA,NA,NA,NA,NA,-2.51,NA),
                            GG=c(NA,NA,NA,NA,NA,-3.26,NA,-1.53,NA,NA,NA,NA,NA,-2.11,NA,-0.5),
                            GU=c(NA,NA,NA,NA,-2.24,NA,-2.51,NA,NA,NA,NA,NA,-1.36,NA,1.29,NA),
                            UA=c(NA,NA,NA,-1.33,NA,NA,NA,NA,NA,NA,NA,-1,NA,NA,NA,NA),
                            UC=c(NA,NA,-2.35,NA,NA,NA,NA,NA,NA,NA,-1.53,NA,NA,NA,NA,NA),
                            UG=c(NA,-2.11,NA,-1,NA,NA,NA,NA,NA,-1.41,NA,0.3,NA,NA,NA,NA),
                            UU=c(-0.93,NA,-1.27,NA,NA,NA,NA,NA,-0.55,NA,-0.5,NA,NA,NA,NA,NA))
  rownames(energyTable) <- colnames(energyTable)
  
  deltaET <- data.frame(AA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.03),
                        AC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.06,NA),
                        AG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.06,NA,0.32),
                        AU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.08,NA,0.24,NA),
                        CA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.07,NA,NA,NA,NA),
                        CC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.07,NA,NA,NA,NA,NA),
                        CG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,0.09,NA,0.24,NA,NA,NA,NA),
                        CU=c(NA,NA,NA,NA,NA,NA,NA,NA,0.06,NA,0.25,NA,NA,NA,NA,NA),
                        GA=c(NA,NA,NA,NA,NA,NA,NA,0.06,NA,NA,NA,NA,NA,NA,NA,0.28),
                        GC=c(NA,NA,NA,NA,NA,NA,0.08,NA,NA,NA,NA,NA,NA,NA,0.25,NA),
                        GG=c(NA,NA,NA,NA,NA,0.07,NA,0.27,NA,NA,NA,NA,NA,0.25,NA,0.96),
                        GU=c(NA,NA,NA,NA,0.06,NA,0.25,NA,NA,NA,NA,NA,0.24,NA,0.56,NA),
                        UA=c(NA,NA,NA,0.09,NA,NA,NA,NA,NA,NA,NA,0.3,NA,NA,NA,NA),
                        UC=c(NA,NA,0.06,NA,NA,NA,NA,NA,NA,NA,0.27,NA,NA,NA,NA,NA),
                        UG=c(NA,0.07,NA,0.3,NA,NA,NA,NA,NA,0.24,NA,0.48,NA,NA,NA,NA),
                        UU=c(0.03,NA,0.28,NA,NA,NA,NA,NA,0.32,NA,0.96,NA,NA,NA,NA,NA))
  rownames(deltaET) <- colnames(deltaET)
  
  endMM <- data.frame(A=c(T,T,T,F), C=c(T,T,F,T), G=c(T,F,T,F), U=c(F,T,F,T))
  rownames(endMM) <- colnames(endMM)
  
  #terminal mismatches Turner 2004 https://rna.urmc.rochester.edu/NNDB/turner04/index.html
  mmNrg <- data.frame(AA=c(NA,NA,NA,-1.0,NA,NA,NA,-0.7,NA,NA,NA,-1.1,-0.8,-1.0,-0.8,-1.0),
                      AC=c(NA,NA,-1.1,NA,NA,NA,-1.1,NA,NA,NA,-1.6,NA,-0.6,-0.7,-0.6,-0.7),
                      AG=c(NA,-1.5,NA,-1.0,NA,-1.0,NA,-0.7,NA,-1.4,NA,-0.5,-0.8,-1.0,-0.8,-1.0),
                      AU=c(-0.8,NA,-0.3,NA,-0.6,NA,-0.6,NA,-0.8,NA,-0.6,NA,-0.6,-0.8,-0.6,-0.8),
                      CA=c(NA,NA,NA,-0.8,NA,NA,NA,-0.6,-1.5,-1.5,-1.4,-1.5,NA,NA,NA,-0.6),
                      CC=c(NA,NA,-1.5,NA,NA,NA,-0.7,NA,-1.0,-1.1,-1.0,-0.8,NA,NA,-1.0,NA),
                      CG=c(NA,-1.5,NA,-0.8,NA,-1.1,NA,-0.6,-1.4,-1.5,-1.6,-1.5,NA,-1.4,NA,-0.6),
                      CU=c(-1.0,NA,-1.0,NA,-0.7,NA,-0.7,NA,-1.0,-1.4,-1.0,-1.2,-0.7,NA,-0.8,NA),
                      GA=c(NA,NA,NA,-1.1,-1.1,-1.5,-1.3,-1.5,NA,NA,NA,-1.2,-0.3,-1.0,-0.8,-1.0),
                      GC=c(NA,NA,-1.3,NA,-1.1,-0.7,-1.1,-0.5,NA,NA,-1.4,NA,-0.6,-0.7,-0.6,-0.7),
                      GG=c(NA,-1.4,NA,-1.1,-1.6,-1.5,-1.4,-1.5,NA,-1.6,NA,-0.8,-0.6,-1.0,-0.8,-1.0),
                      GU=c(-0.8,NA,-0.8,NA,-1.1,-1.0,-1.1,-0.7,-0.8,NA,-0.8,NA,-0.6,-0.8,-0.6,-0.6),
                      UA=c(-1.0,-0.8,-1.1,-0.8,NA,NA,NA,-0.5,-1.0,-0.8,-1.1,-0.8,NA,NA,NA,-0.5),
                      UC=c(-0.7,-0.6,-0.7,-0.5,NA,NA,-0.5,NA,-0.7,-0.6,-0.7,-0.5,NA,NA,-0.7,NA),
                      UG=c(-1.1,-0.8,-1.2,-0.8,NA,-0.8,NA,NA,-0.5,-0.8,-0.8,-0.8,NA,-1.2,NA,-0.5),
                      UU=c(-0.7,-0.6,-0.7,-0.5,-0.7,NA,NA,NA,-0.7,-0.6,-0.7,-0.5,-0.8,NA,-0.6,NA),
                      stringsAsFactors=F)
  rownames(mmNrg) <- colnames(mmNrg)
  
  #dangling ends Turner 2004 https://rna.urmc.rochester.edu/NNDB/turner04/index.html
  dangleStart <- data.frame(A=c(NA,NA,NA,NA,NA,NA,NA,-0.7,NA,NA,NA,-0.1,NA,NA,NA,-0.7,NA,NA,NA,-0.1),	
                            C=c(NA,NA,NA,NA,NA,NA,-1.1,NA,NA,NA,-0.4,NA,NA,NA,-1.3,NA,NA,NA,-0.6,NA),	
                            G=c(NA,NA,NA,NA,NA,-1.7,NA,-0.7,NA,-0.8,NA,-0.1,NA,-1.7,NA,-0.7,NA,-1.2,NA,-0.1),	
                            U=c(NA,NA,NA,NA,-0.8,NA,-0.8,NA,-0.5,NA,-0.5,NA,-0.8,NA,-0.8,NA,-0.6,NA,-0.6,NA),	
                            AA=c(NA,NA,NA,-0.3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AC=c(NA,NA,-0.5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AG=c(NA,-0.2,NA,-0.3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AU=c(-0.3,NA,-0.3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CA=c(NA,NA,NA,-0.3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CC=c(NA,NA,-0.3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CG=c(NA,-0.3,NA,-0.3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CU=c(-0.1,NA,-0.1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GA=c(NA,NA,NA,-0.4,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GC=c(NA,NA,-0.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GG=c(NA,0,NA,-0.4,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GU=c(-0.2,NA,-0.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UA=c(NA,NA,NA,-0.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UC=c(NA,NA,-0.1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UG=c(NA,0,NA,-0.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UU=c(-0.2,NA,-0.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA), 
                            stringsAsFactors=F)
  rownames(dangleStart) <- colnames(dangleStart)
  
  dangleEnds <- data.frame(A=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.3,-0.1,-0.2,-0.2),	
                           C=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.2,-0.3,0,0,NA,NA,NA,NA),	
                           G=c(NA,NA,NA,NA,NA,NA,NA,NA,-0.5,-0.3,-0.2,-0.1,NA,NA,NA,NA,-0.3,-0.1,-0.2,-0.2),	
                           U=c(NA,NA,NA,NA,-0.3,-0.3,-0.4,-0.2,NA,NA,NA,NA,-0.3,-0.3,-0.4,-0.2,NA,NA,NA,NA),	
                           AA=c(NA,NA,NA,-0.8,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           AC=c(NA,NA,NA,-0.5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           AG=c(NA,NA,NA,-0.8,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           AU=c(NA,NA,NA,-0.6,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           CA=c(NA,NA,-1.7,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           CC=c(NA,NA,-0.8,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           CG=c(NA,NA,-1.7,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           CU=c(NA,NA,-1.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           GA=c(NA,-1.1,NA,-0.8,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           GC=c(NA,-0.4,NA,-0.5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           GG=c(NA,-1.3,NA,-0.8,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           GU=c(NA,-0.6,NA,-0.6,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           UA=c(-0.7,NA,-0.7,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           UC=c(-0.1,NA,-0.1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           UG=c(-0.7,NA,-0.7,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                           UU=c(-0.1,NA,-0.1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA), 
                           stringsAsFactors=F)
  rownames(dangleEnds) <- colnames(dangleEnds)
  
  #single nucleotide mismatches Davis 2010 DOI: 10.1021/bi100146z
  mismatchNrg <- data.frame(UAC=c(NA,-0.64,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.77,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AGG=c(-0.64,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.01,NA,NA,NA),	
                            CAU=c(NA,NA,NA,-0.64,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.77,NA,NA,NA,NA,NA,NA,NA),	
                            GGA=c(NA,NA,-0.64,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.01,NA,NA),	
                            CAC=c(NA,NA,NA,NA,NA,0.21,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GGG=c(NA,NA,NA,NA,0.21,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UAG=c(NA,NA,NA,NA,NA,NA,NA,1.26,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AGC=c(NA,NA,NA,NA,NA,NA,1.26,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GAU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,1.26,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CGA=c(NA,NA,NA,NA,NA,NA,NA,NA,1.26,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UAU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.92,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AGA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.92,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UUG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.82,NA,NA,0.33,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GUC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.82,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GUU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.82,NA,0.33,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CUG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-2.82,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AUC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.33,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CUA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.33,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.17,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CAA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.17,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UCG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.17,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AAC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.17,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CAG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.32,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GCC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.32,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GAC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.32,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CCG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.32,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            AAG=c(1.77,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.51,NA,NA,NA,NA,NA,NA,NA,NA),	
                            UCC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.51,NA,NA,NA,NA,NA,NA,NA,2.24,NA),	
                            GAA=c(NA,NA,1.77,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.51,NA,NA,NA,NA,NA,NA),	
                            CCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.51,NA,NA,NA,NA,NA,NA,2.24),	
                            UCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.26,NA,NA,NA,NA),	
                            AUA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.26,NA,NA,NA,NA,NA),	
                            UGC=c(NA,-0.01,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            CGU=c(NA,NA,NA,-0.01,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                            ACG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.24,NA,NA,NA,NA,NA,NA,NA,NA),	
                            GCA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.24,NA,NA,NA,NA,NA,NA), 
                            stringsAsFactors=F)
  rownames(mismatchNrg) <- colnames(mismatchNrg)
  
  #5' shift penalties DOI: 10.1021/bi100146z
  mismatch5Shift <- data.frame(UAC=c(NA,0.425,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.08,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AGG=c(0.24,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.065,NA,NA,NA),	
                               CAU=c(NA,NA,NA,0.24,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.52,NA,NA,NA,NA,NA,NA,NA),	
                               GGA=c(NA,NA,0.425,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.14,NA,NA),	
                               CAC=c(NA,NA,NA,NA,NA,-0.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GGG=c(NA,NA,NA,NA,-0.53,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UAG=c(NA,NA,NA,NA,NA,NA,NA,-1.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AGC=c(NA,NA,NA,NA,NA,NA,-0.18,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GAU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.18,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CGA=c(NA,NA,NA,NA,NA,NA,NA,NA,-1.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UAU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.265,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AGA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.715,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UUG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.63,NA,NA,-0.065,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GUC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.165,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GUU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.165,NA,-0.49,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CUG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.63,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AUC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.49,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CUA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.065,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.265,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CAA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.165,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UCG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.165,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AAC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.265,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CAG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GCC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.185,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GAC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.185,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CCG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AAG=c(-0.52,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.43,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UCC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA,NA,NA,NA,0,NA),	
                               GAA=c(NA,NA,-0.08,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA,NA,NA),	
                               CCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.43,NA,NA,NA,NA,NA,NA,-0.88),	
                               UCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.08,NA,NA,NA,NA),	
                               AUA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.52,NA,NA,NA,NA,NA),	
                               UGC=c(NA,0.14,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CGU=c(NA,NA,NA,-0.065,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               ACG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.88,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GCA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA,NA,NA),
                               stringsAsFactors=F)
  rownames(mismatch5Shift) <- colnames(mismatch5Shift)
  
  #3' shift penalties DOI: 10.1021/bi100146z
  mismatch3Shift <- data.frame(UAC=c(NA,0.24,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.52,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AGG=c(0.425,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.14,NA,NA,NA),	
                               CAU=c(NA,NA,NA,0.425,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.08,NA,NA,NA,NA,NA,NA,NA),	
                               GGA=c(NA,NA,0.24,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.065,NA,NA),	
                               CAC=c(NA,NA,NA,NA,NA,-0.53,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GGG=c(NA,NA,NA,NA,-0.2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UAG=c(NA,NA,NA,NA,NA,NA,NA,-0.18,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AGC=c(NA,NA,NA,NA,NA,NA,-1.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GAU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,-1.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CGA=c(NA,NA,NA,NA,NA,NA,NA,NA,-0.18,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UAU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.715,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AGA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.265,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UUG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.165,NA,NA,-0.49,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GUC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.63,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GUU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.63,NA,-0.065,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CUG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.165,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AUC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.065,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CUA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.49,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.165,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CAA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.265,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UCG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.265,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AAC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.165,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CAG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.185,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GCC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GAC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.025,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CCG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.185,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               AAG=c(-0.08,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA),	
                               UCC=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.43,NA,NA,NA,NA,NA,NA,NA,-0.88,NA),	
                               GAA=c(NA,NA,-0.52,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.43,NA,NA,NA,NA,NA,NA),	
                               CCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA,NA,NA,0),	
                               UCU=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.52,NA,NA,NA,NA),	
                               AUA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.08,NA,NA,NA,NA,NA),	
                               UGC=c(NA,-0.065,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               CGU=c(NA,NA,NA,0.14,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),	
                               ACG=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA),	
                               GCA=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.88,NA,NA,NA,NA,NA,NA),
                               stringsAsFactors=F)
  rownames(mismatch3Shift) <- colnames(mismatch3Shift)
  
  #end penalties Turner 2004 https://rna.urmc.rochester.edu/NNDB/turner04/index.html
  endPenalties <- data.frame(AU=.45, CG=0, GC=0, GU=.45, UA=.45, UG=.45, stringsAsFactors = F)
  deltaEP <- data.frame(AU=.04, CG=0, GC=0, GU=.04, UA=.04, UG=.04, stringsAsFactors = F)
  
  if(calc=='min'){
    energyTable <- energyTable - deltaET
    endPenalties <- endPenalties - deltaEP
  } else if(calc=='max') {
    energyTable <- energyTable + deltaET
    endPenalties <- endPenalties + deltaEP
  }
  
  #make data frame
  if(is.na(x53) && is.na(y35)){ #sample data
    df <- data.frame(seq=c("G", 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA','GCAAAAG', 'GCGAAAG', 'GCAAAAGC', 'GCGAAAGC', "G", 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA','GCAAAAG', 'GCGAAAG', 'GCAAAAGC', 'GCGAAAGC'), template=c('UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC','UCGUUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC', 'UCGCUUUCGUCC'), stringsAsFactors = F)
  } else {
    if(length(x53)!=length(y35)){
      warning("different numbers of sequences in each group")
    }
    df <- data.frame(seq=x53, template=y35, stringsAsFactors = F) #passed data
  }
  
  #offset handling
  if(length(offset)==1){ #works with a vector or an int
    if(offset>0){
      
      adj <- " "
      if(offset>1){
        for(i in 2:offset){
          adj <- paste(adj," ", sep="")
        }
      }
      
      df$seq <- paste(adj, df$seq, sep="")
      
    } else if(offset<0){
      offset <- offset*-1
      
      adj <- " "
      if(offset>1){
        for(i in 2:offset){
          adj <- paste(adj," ", sep="")
        }
      }
      
      #offset testing
      #  nums <- 1:10
      #  adj <- nums[1]
      #  if(offset>1){
      #    for(i in 2:offset){
      #      adj <- paste(adj,nums[i], sep="")
      #    }
      #  }
      
      df$template <- paste(adj, df$template, sep="")
    }
    
    offset <- offset + 1 #fixes offset to new start position
    df$startPos <- offset
    
  } else if(length(offset)==nrow(df)){
    df$startPos <- offset
    df$offset <- ""
    for(i in 1:max(abs(df$startPos))){
      df$offset[abs(df$startPos)>=i] <- paste(df$offset[abs(df$startPos)>=i], " ", sep="")
    }
    df$seq[df$startPos>0] <- paste(df$offset[df$startPos>0], df$seq[df$startPos>0], sep="")
    df$template[df$startPos<0] <- paste(df$offset[df$startPos<0], df$template[df$startPos<0], sep="")
    df$startPos[df$startPos<0] <- df$startPos[df$startPos<0] *-1
    
    df <- df[,-ncol(df)]
    df$startPos <- df$startPos + 1
  } else {
    stop("Offset length must be a single int or a vector with as many elements as x53 and y35")
  }
  
  #nchar stop point for each row
  pnums <- data.frame(nchar(df[,1]),nchar(df[,2]))
  pnums[,3] <- pnums[,1]
  pnums[,3][pnums[,2]<pnums[,1]] <- pnums[,2][pnums[,2]<pnums[,1]]
  df$endPos <- as.vector(pnums[,3])
  rm(pnums)
  
  #start and end of helix
  df$start <- paste(substr(df$seq, df$startPos, df$startPos),substr(df$template, df$startPos, df$startPos), sep="-")
  df$end <- paste(substr(df$seq, df$endPos, df$endPos),substr(df$template,df$endPos, df$endPos), sep="-")
  
  #dangling ends
  df$dangle53 <- paste(substr(df$seq, df$startPos-1, df$startPos),substr(df$template, df$startPos-1, df$startPos), sep="-")
  df$dangle53 <- gsub(" ","",df$dangle53)
  df$dangle53[nchar(df$dangle53)!=4] <- ""
  
  df$dangle35 <- ""
  df$dangle35[nchar(df$seq)>nchar(df$template)] <- paste(substr(df$seq[nchar(df$seq)>nchar(df$template)], nchar(df$template[nchar(df$seq)>nchar(df$template)]), nchar(df$template[nchar(df$seq)>nchar(df$template)])+1),substr(df$template[nchar(df$seq)>nchar(df$template)], nchar(df$template[nchar(df$seq)>nchar(df$template)]), nchar(df$template[nchar(df$seq)>nchar(df$template)])), sep="-")
  df$dangle35[nchar(df$seq)<nchar(df$template)] <- paste(substr(df$seq[nchar(df$seq)<nchar(df$template)], nchar(df$seq[nchar(df$seq)<nchar(df$template)]), nchar(df$seq[nchar(df$seq)<nchar(df$template)])),substr(df$template[nchar(df$seq)<nchar(df$template)], nchar(df$seq[nchar(df$seq)<nchar(df$template)]), nchar(df$seq[nchar(df$seq)<nchar(df$template)])+1), sep="-")
  df$dangle35 <- gsub(" ","",df$dangle35)
  df$dangle35[nchar(df$dangle35)!=4] <- ""
  
  #finding mismatches at either end
  df$startMM <- ""
  df$endMM <- ""
  repeat{
    df$sMM <- F
    df$eMM <- F
    
    for(i in 1:ncol(endMM)){
      for(j in 1:nrow(endMM)){
        end <- paste(colnames(endMM)[i], rownames(endMM)[j], sep="-")
        df$sMM[df$start==end] <- endMM[j,i]
        df$eMM[df$end==end] <- endMM[j,i]
      }
      rm(end)
    }
    
    if(all(df$sMM==F) && all(df$eMM==F)){
      break
    }
    
    #if mismatch, remove dangling ends, get mismatch, set new start/end position, change start and end sequence
    df$dangle53[df$sMM] <- ""
    df$startMM[df$sMM] <- paste(substr(df$seq[df$sMM], df$startPos[df$sMM], df$startPos[df$sMM]+1),substr(df$template[df$sMM], df$startPos[df$sMM], df$startPos[df$sMM]+1), sep="-")
    df$startPos[df$sMM] <- df$startPos[df$sMM] + 1
    
    df$dangle35[df$eMM] <- ""
    df$endMM[df$eMM] <- paste(substr(df$seq[df$eMM], df$endPos[df$eMM]-1, df$endPos[df$eMM]),substr(df$template[df$eMM], df$endPos[df$eMM]-1, df$endPos[df$eMM]), sep="-")
    df$endPos[df$eMM] <- df$endPos[df$eMM] - 1
    
    df$start <- paste(substr(df$seq, df$startPos, df$startPos),substr(df$template, df$startPos, df$startPos), sep="-")
    df$end <- paste(substr(df$seq, df$endPos, df$endPos),substr(df$template,df$endPos, df$endPos), sep="-")
  } #check fdor more mismatches
  
  df <- df[,-c(ncol(df)-1, ncol(df))]
  
  #nearest neighbour pairs
  for(i in 1:(max(df$endPos - df$startPos))){
    df[,ncol(df)+1] <- ""
    df[,ncol(df)][df$endPos>=i+df$startPos] <- paste(substr(df$seq[df$endPos>=i+df$startPos], df$startPos[df$endPos>=i+df$startPos]+i-1, df$startPos[df$endPos>=i+df$startPos]+i), substr(df$template[df$endPos>=i+df$startPos], df$startPos[df$endPos>=i+df$startPos]+i-1, df$startPos[df$endPos>=i+df$startPos]+i), sep="-")
    df[,ncol(df)][df$endPos<i+df$startPos] <- ""
    colnames(df)[ncol(df)] <- paste("pair",i,sep="_")
    pairs <- i
  }
  
  #convert to energy
  for(i in 1:ncol(energyTable)){
    for(j in 1:nrow(energyTable)){
      df[,(ncol(df)-pairs+1):ncol(df)][df[,(ncol(df)-pairs+1):ncol(df)]==paste(colnames(energyTable)[i], rownames(energyTable)[j], sep="-")] <- energyTable[j,i]
    }
  }
  if(pairs>1){
    df[rowSums(is.na(df[,(ncol(df)-pairs+1):ncol(df)]))==1,][is.na(df[rowSums(is.na(df[,(ncol(df)-pairs+1):ncol(df)]))==1,])] <- 0
  }
  df[,c(grep("pair_", colnames(df)))][df[,c(grep("pair_", colnames(df)))]==""] <- 0
  
  
  ####Start NA handling
  nas <- max(rowSums(is.na(df)))
  if(nas>0){ #slower if nas, but faster if no mismatch
    
    #find all the nas, label as T in new columns
    df[,(ncol(df)+1):(ncol(df)+pairs)] <- is.na(df[,(ncol(df)-pairs+1):ncol(df)])
    df[,(ncol(df)-pairs+1):ncol(df)][df[,(ncol(df)-pairs+1):ncol(df)]==F] <- "F"
    df[,(ncol(df)-pairs+1):ncol(df)][df[,(ncol(df)-pairs+1):ncol(df)]==T] <- "T"
    colnames(df)[(ncol(df)-pairs+1):ncol(df)] <- paste("NAnum_",1:pairs,sep="")
    
    #collapse rows to a string
    pasteArgs <- c(df[,c(grep("NAnum_", colnames(df)))], sep="")
    df[,ncol(df)+1] <- do.call(paste, pasteArgs)
    rm(pasteArgs)
    df <- df[-c(grep("NAnum_", colnames(df)))]
    
    colnames(df)[ncol(df)] <- "ntmm"
    
    for(i in 1:ceiling(nas/2)){ #gets mismatch location by regexpr and collapsed rows string
      df[,ncol(df)+1] <- regexpr("T", df$ntmm)
      colnames(df)[ncol(df)] <- paste("mm",i,sep="_")
      df$ntmm[df[,ncol(df)]>0] <- substr(df$ntmm[df[,ncol(df)]>0], df[,ncol(df)][df[,ncol(df)]>0]+2, nchar(df$ntmm[df[,ncol(df)]>0]))
      df[,ncol(df)][df[,ncol(df)]>0] <- df[,ncol(df)][df[,ncol(df)]>0] + df$startPos[df[,ncol(df)]>0]
    }
    
    #fix up the pair energies
    df[,c(grep("mm_", colnames(df)))][df[,c(grep("mm_", colnames(df)))]<0] <- 0
    df[,c(grep("pair_", colnames(df)))][is.na(df[,c(grep("pair_", colnames(df)))])] <- 0
    df <- df[,-c(grep("ntmm", colnames(df)))]
    
    #get mm substrs
    for(i in 1:ceiling(nas/2)){
      df[,ncol(df)+1] <- ""
      #mm_ gives na location; -1 +1 gives 3 nt needed for snmm
      df[,ncol(df)][df[,grep("mm_",colnames(df))[i]]>0] <- paste(substr(df$seq[df[,grep("mm_",colnames(df))[i]]>0], df[,grep("mm_",colnames(df))[i]][df[,grep("mm_",colnames(df))[i]]>0]-1, df[,grep("mm_",colnames(df))[i]][df[,grep("mm_",colnames(df))[i]]>0]+1), substr(df$template[df[,grep("mm_",colnames(df))[i]]>0], df[,grep("mm_",colnames(df))[i]][df[,grep("mm_",colnames(df))[i]]>0]-1, df[,grep("mm_",colnames(df))[i]][df[,grep("mm_",colnames(df))[i]]>0]+1), sep="-")
      colnames(df)[ncol(df)] <- paste('snMis', i, sep="_")
      
      df[,ncol(df)+1] <- 0
      df[,ncol(df)][df[,grep("mm_",colnames(df))[i]]>0] <- 4 - (df[,grep("mm_",colnames(df))[i]][df[,grep("mm_",colnames(df))[i]]>0] - df$startPos[df[,grep("mm_",colnames(df))[i]]>0])
      df[,ncol(df)][ df[,ncol(df)]<0] <- 0
      colnames(df)[ncol(df)] <- paste('snMisls', i, sep="_")
      
      df[,ncol(df)+1] <- 0
      df[,ncol(df)][df[,grep("mm_",colnames(df))[i]]>0] <- 4 - (df$endPos[df[,grep("mm_",colnames(df))[i]]>0] - df[,grep("mm_",colnames(df))[i]][df[,grep("mm_",colnames(df))[i]]>0])
      df[,ncol(df)][ df[,ncol(df)]<0] <- 0
      colnames(df)[ncol(df)] <- paste('snMisrs', i, sep="_")
    }
    
    #for in nas/2
    #get match
    #check 5' and adj
    #check 3' and adj
    #get the energy of the mismatch
    for(k in 1:ceiling(nas/2)){
      for(i in 1:ncol(mismatchNrg)){
        for(j in 1:nrow(mismatchNrg)){
          df[,grep("snMis_", colnames(df))[k]][df[,grep("snMis_", colnames(df))[k]]==paste(colnames(mismatchNrg)[i],rownames(mismatchNrg)[j], sep="-")] <- mismatchNrg[j,i] +  df[,grep("snMisls_", colnames(df))[k]][df[,grep("snMis_", colnames(df))[k]]==paste(colnames(mismatchNrg)[i],rownames(mismatchNrg)[j], sep="-")] * mismatch5Shift[j,i] + df[,grep("snMisrs_", colnames(df))[k]][df[,grep("snMis_", colnames(df))[k]]==paste(colnames(mismatchNrg)[i],rownames(mismatchNrg)[j], sep="-")] * mismatch3Shift[j,i] 
        }
      }
    }
    df <-  df[,-c(grep("mm_|snMisls|snMisrs",colnames(df)))] #get rid of all added columns
  }#end NA handling
  
  ####start and end handling
  df$start <- gsub("-","", df$start)
  df$end <- gsub("-","", df$end)
  for(i in 1:ncol(endPenalties)){
    df$start[df$start==colnames(endPenalties)[i]] <- endPenalties[1,i]
    df$end[df$end==colnames(endPenalties)[i]] <- endPenalties[1,i]
  }
  
  ####dangling ends
  for(i in 1:ncol(dangleStart)){
    for(j in 1:nrow(dangleStart)){
      df$dangle53[df$dangle53==paste(colnames(dangleStart)[i],rownames(dangleStart)[j],sep='-')] <- dangleStart[j,i]
      df$dangle35[df$dangle35==paste(colnames(dangleEnds)[i],rownames(dangleEnds)[j],sep='-')] <- dangleEnds[j,i]
    }
  }
  
  ####terminal mismatches
  for(i in 1:ncol(mmNrg)){
    for(j in 1:nrow(mmNrg)){
      df$startMM[df$startMM==paste(colnames(mmNrg)[i], rownames(mmNrg)[j], sep="-")] <- mmNrg[j,i]
      df$endMM[df$endMM==paste(colnames(mmNrg)[i], rownames(mmNrg)[j], sep="-")] <- mmNrg[j,i]
    }
  }
  
  #setting zeros, converting data to numeric, taking sum of data
  df[,5:ncol(df)][df[,5:ncol(df)]==""] <- 0
  df[,5:ncol(df)] <- as.numeric(as.character(unlist(df[,5:ncol(df)])))
  df$totalEnergy <- rowSums(df[,5:ncol(df)])
  
  #df.helnrg <<- df
  
  return(df$totalEnergy)
}
strBraceRemoval <- function(x) {
  # Removes braces from a string and returns a vector of strings containing all possible strings with
  # the contained characters. i.e., 12[34] would return 123 and 124
  #
  # Args:
  #   x: a string or vector of strings
  #
  # Returns:
  #   strBraceRemoval(str): recursion, removes subsequent []
  #   str: vector containing strings with [] removed and replaced with contents 
  
  #take string
  #find [
  #extract things between
  #get number
  #recompile string with each possibility
  #if it finds more []
  #recursion
  x <- as.vector(x)
  str <- vector(mode="character", length=0)
  for(i in 1:length(x)){
    if(grepl("\\[", x[i])){
      x.open <- regexpr("\\[", x[i], perl=T)
      x.close <- regexpr("\\]", x[i], perl=T)
      x.contents <- substr(x[i], x.open+1, x.close-1)
      x.contents <- unlist(strsplit(x.contents, split=""))
      for(j in 1:length(x.contents)){
        str[(length(str)+1)] <- paste(substr(x[i], 1, x.open-1),x.contents[j], substr(x[i], x.close+1, nchar(x[i])), sep="")
      }
    } else {
      str[(length(str)+1)] <- x[i]
    }
  }
  if(any(grepl("\\[",str))){
    return(strBraceRemoval(str))
  } else {
    return(str)
  }
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
interLength <- function(df, rounds, trim, seqAdded=NA, seqRound, numRound=NA, final=F){
  #takes a dataframe of RNA reads that have been trimmed for prime and realign, gives intermediates between min rounds + 1 and max rounds - 1
  #
  #Args:
  # df: dataframe containing trimmed sequence, rounds, and additions
  # rounds: column index containing rounds of trimming information
  # trim: column index containing the trimmed sequence
  # seqAdded: column index containing the sequence added (optional)
  # seqRound: column index containing the seq added each round; may be a vector
  # numRound: column index containing the number of nucleotides added each round; may be a vector
  #   optional argument will also recalculate other numeric markers
  # final: include the fully realigned intermediate
  #
  #return:
  # df: dataframe containing the original with the intermediate trimming products added
  
  if(final==F){
    loopVar <- (max(df[,rounds])-1) #max rounds -1 is the upper limit for intermediates
  } else {
    loopVar <- (max(df[,rounds]))
  }
  df$Intermediate <- FALSE
  
  for(i in 1:loopVar){
    if(i==1){ #creates temporary df, or trims existing temporary df
      dfTemp <- df[df[,rounds]>1,]
      dfTemp$Intermediate <- TRUE
    } else if(final && loopVar==i){
      dfTemp <- df[df[,rounds]==1,]
      dfTemp$Intermediate <- TRUE
    } else {
      dfTemp <- dfTemp[dfTemp[,rounds]>1,]
    }
    
    dfTemp[,rounds] <- dfTemp[,rounds]-1 #drops rounds by one so as to show number of rounds trimming to the new length
    dfTemp[,trim] <- paste(dfTemp[,trim],dfTemp[,seqRound[1]],sep="") #adds the removed sequence
    
    for(j in 1:loopVar){ #left shifts the sequence added
      if(j!=length(seqRound)){
        dfTemp[,seqRound[j]] <- dfTemp[,seqRound[j+1]]
        dfTemp[,seqRound[j+1]] <- ""
      } else {
        dfTemp[,seqRound[j]] <- ""
      }
    }
    
    if(!is.na(seqAdded)){ #changes what was added to reflect intermediate
      paste_args <- c(dfTemp[,c(grep("Seq_Trim_", colnames(dfTemp)))], sep="")
      dfTemp[,seqAdded] <- do.call(paste,paste_args)
    }
    
    df <- rbind(df, dfTemp) #adds intermediates
  }
  
  if(!is.na(numRound[1])){ #checks to see if numbers should be addressed by this program
    for(i in 1:length(numRound)){ #fixes numbers
      df[,numRound[i]] <- nchar(df[,seqRound[i]])
    }
    
    #adds/changes number variables
    df$Trim_Len <- nchar(df[,trim])
    df$Len_Adjust <- nchar(df[,seqAdded])
    df$Velocity <- df$Len_Adjust/df[,rounds]
    df$Velocity[is.na(df$Velocity)] <- 0
  }
  dfintlenfun <<- df
  return(df)
}

#setwd("") #if not default dir
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
wsn <- read.csv("WSN_All_Match.csv", stringsAsFactors = F)
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)

ntCap <- 3

####Generate conversion table
templates <- c("UCGUUUUCGUCC", "UCGCUUUCGUCC")
additionVect <- c("[ACGU]", '[ACGU]', 'G','C','[AG]','A','A','A','G','C','A','G','G') #string objects containing all possibilities
additions <- strBraceRemoval(paste(additionVect[1:3],collapse = ""))
if(length(additionVect)>=4){
  for(i in 4:length(additionVect)){
    additions <- c(additions, strBraceRemoval(paste(additionVect[1:i],collapse = ""))) #increase vector length 
  }
}
rm(additionVect)

i=1
if(length(additions) <= length(templates)){
  dfArgs <- list(Additions=additions[i], Templates=templates, stringsAsFactors=F)
} else {
  dfArgs <- list(Additions=additions, Templates=templates[i], stringsAsFactors=F)
}

nrgConverter <- do.call(data.frame, dfArgs)
rm(dfArgs)

if(min(length(additions), length(templates))>1){
  for(i in 2:min(length(additions), length(templates))){
    if(length(additions) <= length(templates)){
      dfArgs <- list(Additions=additions[i], Templates=templates, stringsAsFactors=F)
    } else {
      dfArgs <- list(Additions=additions, Templates=templates[i], stringsAsFactors=F)
    }
    nrgConverter <- rbind(nrgConverter, do.call(data.frame, dfArgs))
    rm(dfArgs)
  }
}
rm(i, templates, additions)
nrgConverter$TotalEnergy <- helixEnergy(x53=nrgConverter[,1], y35=nrgConverter[,2], offset=-1, calc='max') 
####End conversion table generation

dfPreProcess <- function(df, subsa=c(1:5), subsb=c(6:8), converter=NA, strain="Influenza"){
  #pre-processes data for bon energy.R
  #
  #args
  # df: dataframe
  # subs: columns containing subunit reads
  # converter:  an nrgConverter made by helixEnergy for all possible templates and additions; if NA returns the table
  # strain: name of the strain; stored in a column for global analysis
  #
  #return
  # df: processed data.frame
  
  subs <- c(subsa, subsb)
  
  df <- df[,c(subs, grep("round", colnames(df)), grep("Trim_Sequence",colnames(df)), grep("Seq_Trim_R", colnames(df)), grep("NT_coupure", colnames(df)), grep("Sequence",colnames(df))[1])]
  df <- interLength(df, (length(subs)+1),grep("Trim_Sequence",colnames(df)), seqRound=grep("Seq_Trim_R", colnames(df)), final = T)
  df$Strain <- strain
  df <- df[,c(subs, ncol(df), grep("Trim_Sequence",colnames(df)), grep("Seq_Trim_R", colnames(df))[1], grep("NT_coupure", colnames(df)),(ncol(df)-1))]
  
  df$NT_coupure <- substr(df$NT_coupure, 1, 1)
  df <- subunitSeperate(df, subs)
  
  df$Template <- ""
  df$Template[rowSums(df[,c(subsa)])>0] <- "UCGUUUUCGUCC"
  df$Template[rowSums(df[,c(subsb)])>0] <- "UCGCUUUCGUCC"
  
  df$Addition <- ""
  df$Addition[df$Seq_Trim_R1!=""] <- paste(substr(df$Trim_Sequence[df$Seq_Trim_R1!=""], nchar(df$Trim_Sequence[df$Seq_Trim_R1!=""])-1, nchar(df$Trim_Sequence[df$Seq_Trim_R1!=""])), df$Seq_Trim_R1[df$Seq_Trim_R1!=""], sep="")
  df$Addition[df$Addition=="" & df$Template=="UCGUUUUCGUCC"] <- paste(substr(df$Trim_Sequence[df$Addition=="" & df$Template=="UCGUUUUCGUCC"], nchar(df$Trim_Sequence[df$Addition=="" & df$Template=="UCGUUUUCGUCC"])-1, nchar(df$Trim_Sequence[df$Addition=="" & df$Template=="UCGUUUUCGUCC"])), "GCAAAAGCAGG", sep="")
  df$Addition[df$Addition=="" & df$Template=="UCGCUUUCGUCC"] <- paste(substr(df$Trim_Sequence[df$Addition=="" & df$Template=="UCGCUUUCGUCC"], nchar(df$Trim_Sequence[df$Addition=="" & df$Template=="UCGCUUUCGUCC"])-1, nchar(df$Trim_Sequence[df$Addition=="" & df$Template=="UCGCUUUCGUCC"])), "GCGAAAGCAGG", sep="")
  df$Addition <- gsub("T","U",df$Addition)
  df$G_comp <- F
  df$G_comp[df$NT_coupure=="G" & df$Intermediate==F] <- T
  df$Realign <- nchar(df$Seq_Trim_R1)
  df$Length <- nchar(df$Trim_Sequence)
  df <- df[,-c(grep("NT_coupure", colnames(df)), grep("Seq_Trim_R", colnames(df)), grep("Trim_Sequence",colnames(df)))]
  
  if(all(is.na(converter))){
    return(df)
  }
  
  df$TotalEnergy <- 0
  
  for(j in 1:(max(nchar(converter$Additions))-2)){
    df[,ncol(df)+1] <- 0
    for(i in 1:nrow(converter[nchar(converter$Additions)==(j+2),])){
      df[,ncol(df)][substr(df$Addition,1,nchar(converter[nchar(converter$Additions)==(j+2),][i,1]))==converter[nchar(converter$Additions)==(j+2),][i,1] & df$Template==converter[nchar(converter$Additions)==(j+2),][i,2]] <- converter[nchar(converter$Additions)==(j+2),][i,3]
      df$TotalEnergy[df$Addition==converter[nchar(converter$Additions)==(j+2),][i,1] & df$Template==converter[nchar(converter$Additions)==(j+2),][i,2]] <- converter[nchar(converter$Additions)==(j+2),][i,3]
    }
    colnames(df)[ncol(df)] <- paste("nt_addition_",j,sep="")
  }
  
  df <- subunitSeperate(df, subs=subs, decompress = T)
  #df[,9:20] <- -1*log2((-1*df[9:20]))
  #df[,9:20][df[,9:20]==Inf] <- 0
  
  return(df)
}
#5x5 graph looks good with 1.6 for legend
#may want to show realignment thresholds with dashed lines and the difference between this line and those that go on to the next length
nrgGraph <- function(df=NA, ntCap=NA, name=""){
  #makes a boxplot showing energy values for passed data.frame pre-processed by dfPreProcess in this script
  #
  #args
  # df: a data.frame from dfPreProcess; if NA spoofs the legend for cowplot
  # ntCap:  upper limitt for nt addition to be graphed
  #
  #returns
  # get_legend(nrgGraph): spoofed legend for cowplot
  # NULL: the df had no rows; returns a NULL so it doesn't break the code with an error
  # nrgGraph: boxplot of energies
  
  if(length(df)==1){ #catches df=NA; returns a spoofed legend
    dfpoint <- data.frame(Series=c("No_Realignment", 'Trimmed','Realigns'), y=-1, x=c(1,2,1), weight=1)
    
    nrgGraph <- ggplot() +
      geom_point(data=dfpoint, aes(x=x, y=y, group=Series, alpha=weight, colour=Series),position = 'jitter', show.legend = T) +
      scale_colour_manual(name='legend', breaks=c('No_Realignment', 'Trimmed', 'Realigns'), 
                          values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Trimmed='#D98000', Realigns=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65)),
                          labels=c(No_Realignment='Transcribed', Trimmed='Will Realign', Realigns='Realigns')) +
      scale_alpha(guide=F) +
      theme(legend.title = element_blank(), legend.text = element_text(size=7))
    
    return(get_legend(nrgGraph))  #end spoofed legend
  } else if(nrow(df)==0) {
    return(NULL)
  }
  
  if(is.na(ntCap)){
    ntCap <- length(grep("nt_addition_", colnames(df))) #sets ntCap to maximum addition
  }
  
  df$Series <- ""
  df$Series[df$Realign==0] <- "No_Realignment"
  df$Series[df$Series==""] <- "Trimmed"
  
  nrgGraph <- ggplot() #shell for graph
  
  i = ntCap # lazy version
  
  
  #I need to get the the data for each unique nrg in each of the 3 groups
  #I need to get the mean and sd
  #I need to recompress the reads for alpha values
  print(paste("nt_addition_", i, sep=""))
  
  #make a data.frame with the means and the standard deviations
  stats <- data.frame(Series=c('No_Realignment', 'Trimmed', 'Realigns'), 
                      mean=c(mean(df[,grep("nt_addition_", colnames(df))[i]][df$Realign==0]), mean(df[,grep("nt_addition_", colnames(df))[i]][df$Realign>i]), mean(df[,grep("nt_addition_", colnames(df))[i]][df$Realign==i])),
                      sd=c(sd(df[,grep("nt_addition_", colnames(df))[i]][df$Realign==0]), sd(df[,grep("nt_addition_", colnames(df))[i]][df$Realign>i]), sd(df[,grep("nt_addition_", colnames(df))[i]][df$Realign==i])),
                      x=c(3*(i-1) + 1, x=3*(i-1) + 2, x=3*(i-1) + 3))
  #need to deal with NaN and NA
  
  #write stats
  write.csv(stats, paste(name, " ", i, ' nt', '.csv', sep=''))
  
  nrgs <- unique(df[,grep("nt_addition_", colnames(df))[i]])
  nrgs <- nrgs[order(nrgs, decreasing = T)]
  
  dfPoints <- rbind(
    data.frame(Series='No_Realignment', y=nrgs, number=0, x=3*(i-1) + 1, stringsAsFactors = F),
    data.frame(Series='Trimmed', y=nrgs, number=0, x=3*(i-1) + 2, stringsAsFactors = F),
    data.frame(Series='Realigns', y=nrgs, number=0, x=3*(i-1) + 3, stringsAsFactors = F)
  )
  
  for(j in 1:length(nrgs)){
    dfPoints[(1+(j-1)),3] <- nrow(df[df$Realign==0 & df[,grep("nt_addition_", colnames(df))[i]]==nrgs[j],])
    dfPoints[(length(nrgs)+1+(j-1)),3] <- nrow(df[df$Realign>i & df[,grep("nt_addition_", colnames(df))[i]]==nrgs[j],])
    dfPoints[(2*length(nrgs)+1+(j-1)),3] <- nrow(df[df$Realign==i & df[,grep("nt_addition_", colnames(df))[i]]==nrgs[j],])
  }
  
  dfPoints <- dfPoints[dfPoints$number>0,]
  
  if(nrow(dfPoints[dfPoints$Series=='No_Realignment',])==0){
    dfPoints[nrow(dfPoints)+1,] <- data.frame(Series='No_Realignment', y=10, number=0, x=-10, stringsAsFactors = F)
  }
  
  if(nrow(dfPoints[dfPoints$Series=='Trimmed',])==0){
    dfPoints[nrow(dfPoints)+1,] <- data.frame(Series='Trimmed', y=10, number=0, x=-10, stringsAsFactors = F)
  }
  
  if(nrow(dfPoints[dfPoints$Series=='Realigns',])==0){
    dfPoints[nrow(dfPoints)+1,] <- data.frame(Series='Realigns', y=10, number=0, x=-10, stringsAsFactors = F)
  }
  
  dfPoints$Series <- factor(dfPoints$Series, levels =c('No_Realignment', 'Trimmed', 'Realigns'))
  dfPoints <- dfPoints[order(dfPoints$Series),]
  dfPoints$number[dfPoints$Series=='No_Realignment'] <- dfPoints$number[dfPoints$Series=='No_Realignment'] /max(dfPoints$number[dfPoints$Series=='No_Realignment'])
  dfPoints$number[dfPoints$Series=='Trimmed'] <- dfPoints$number[dfPoints$Series=='Trimmed'] /max(dfPoints$number[dfPoints$Series=='Trimmed'])
  dfPoints$number[dfPoints$Series=='Realigns'] <- dfPoints$number[dfPoints$Series=='Realigns'] /max(dfPoints$number[dfPoints$Series=='Realigns'])
  
  nrgGraph <-nrgGraph +
    geom_point(data=dfPoints, aes(y=y, x=x, group=Series, colour=Series, alpha=number), show.legend=F, position = position_jitter(width=.2)) +
    geom_errorbar(data=stats, aes(ymin=mean-2*sd, ymax=mean+2*sd, x=x)) +
    stat_summary(fun.y='mean', data=stats, colour="#000000", geom="errorbar", aes(group=Series,x=x, y=mean,ymax=..y.., ymin=..y..), width=.75, linetype="dashed", show.legend=F)
  
  
  rm(nrgs, stats, dfPoints)
  
  
  
  #make the graph look nice
  nrgGraph <- nrgGraph +
    labs(y=bquote(Delta~G[37]^{o}*' (kcal/mol)'), x=paste(ntCap, "Nucleotide Sequence Addition", sep=' ')) + #delta G degree 37 and addition
    scale_x_continuous(breaks=c(((ntCap-1)*3+1):(3*ntCap)), labels=c('Transcribed', 'Will Realign', 'Realigns'), limits=c((ntCap-1)*3+.5, 3*ntCap+.5)) +
    scale_y_continuous(limits = c(-8, -5), breaks=seq(-5,-8, length.out = 7), labels= seq(-5,-8, length.out = 7)) +
    scale_colour_manual(name='legend', breaks=c('No_Realignment', 'Trimmed', 'Realigns'), 
                        values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65), Trimmed='#D98000', Realigns=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65)),
                        labels=c(No_Realignment='Transcribed', Trimmed='Will Realign', Realigns='Realigns'))
    theme(legend.title = element_blank()) +
    theme_bw()
  return(nrgGraph)
}
getNrgGraphs <- function(df, ntCap=NA, target='Subunit', seqLen=NA){
  #generates and saves graphs from a df from dfPreProcess()
  #
  #Args
  # df: dataframe from dfPreProcess
  # ntCap:  upper limitt for nt addition to be graphed
  # target: column to create graphs around (ex., Strain, Subunit)
  # seqLen: vector of sequences to examine; if NA analyzes all lengths
  
  colnames(df)[grep(target, colnames(df))] <- 'target'
  subs <- unique(df$target)
  if(target=="Strain"){
    name <- 'Global'
  } else {
    name <- df$Strain[1]
  }
  
  if(length(seqLen)==1 && is.na(seqLen)){
    seqLen <- c(min(df$Length):max(df$Length))
  }
  
  if(!is.na(ntCap)){
    df <- df[,-c(grep("nt_addition_", colnames(df))[-c(1:ntCap)])]
  } else {
    ntCap <- length(grep("nt_addition_", colnames(df)))
  }
  
  #LoN <- nrgGraph() #Legend of Energy: the Graph Legend
  
  print(target)
  
  lenGraph <- function(df, var, ntCap, sub, target, name){
    LoN <- nrgGraph() #Legend of Energy: the Graph Legend
    for(j in 1:ceiling(length(var)/4)){
      
      
      print(paste("Length =",var[4*(j-1)+1]))
      titleText <- gsub("\\.","",paste(sub," ", target, " - ", var[4*(j-1)+1], " Nucleotides", sep=""))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph1 <- plot_grid(title, nrgGraph(df[df$Length==var[4*(j-1)+1],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
      
      if(length(var)<(4*(j-1)+2)){
        graph2 <- NULL
      } else {
        print(paste('length =',var[4*(j-1)+2]))
        titleText <- gsub("\\.","",paste(sub," ", target, " - ", var[4*(j-1)+2], " Nucleotides", sep=""))
        title <- ggdraw() + draw_label(titleText, fontface='bold')
        graph2 <- plot_grid(title, nrgGraph(df[df$Length==var[4*(j-1)+2],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
      }
      
      if(length(var)<(4*(j-1)+3)){
        graph3 <- NULL
      } else {
        print(paste('length =',var[4*(j-1)+3]))
        titleText <- gsub("\\.","",paste(sub," ", target, " - ", var[4*(j-1)+3], " Nucleotides", sep=""))
        title <- ggdraw() + draw_label(titleText, fontface='bold')
        graph3 <- plot_grid(title, nrgGraph(df[df$Length==var[4*(j-1)+3],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
      }
      
      if(length(var)<(4*(j-1)+4)){
        graph4 <- NULL
      } else {
        print(paste('length =',var[4*(j-1)+4]))
        titleText <- gsub("\\.","",paste(sub," ", target, " - ", var[4*(j-1)+4], " Nucleotides", sep=""))
        title <- ggdraw() + draw_label(titleText, fontface='bold')
        graph4 <- plot_grid(title, nrgGraph(df[df$Length==var[4*(j-1)+4],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
      }
      
      if(j==1){
        graph <- plot_grid(graph1, graph2, graph3, graph4, ncol=4, nrow=1, rel_widths = c(5,5,5,5), labels = c(LETTERS[(4*(j-1)+1):(4*(j-1)+4)], NULL))
      } else {
        graph <- plot_grid(graph, plot_grid(graph1, graph2, graph3, graph4, ncol=4, nrow=1, rel_widths = c(5,5,5,5), labels = c(LETTERS[(4*(j-1)+1):(4*(j-1)+4)], NULL)), rel_widths=c(1,1), rel_heights=c(j-1, 1), ncol=1, nrow=2)
      }
      
      rm(graph1, graph2, graph3, graph4)
      gc()
      
      if(j==ceiling(length(var)/4)){
        print('plotting')
        save_plot(paste(gsub("\\.","",paste(name,"_", sub, "_by_Length_",sep="")),'.png',sep=""), graph, base_height = 5.5*j, base_width = 20, dpi=600)
      }
    }
  }
  
  for(i in 1:ceiling(length(subs)/4)){
    print(subs[4*(i-1)+1])
    lenGraph(df[df$target==subs[4*(i-1)+1],], seqLen, ntCap, sub=subs[4*(i-1)+1], name=name, target=target)
    titleText <- gsub("\\.","",paste(subs[4*(i-1)+1],target, sep=" "))
    title <- ggdraw() + draw_label(titleText, fontface='bold')
    graph1 <- plot_grid(title, nrgGraph(df[df$target==subs[4*(i-1)+1],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
    
    if(length(subs)<(4*(i-1)+2)){
      graph2 <- NULL
    } else {
      print(subs[4*(i-1)+2])
      lenGraph(df[df$target==subs[4*(i-1)+2],], seqLen, ntCap, sub=subs[4*(i-1)+2], name=name, target=target)
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+2],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph2 <- plot_grid(title, nrgGraph(df[df$target==subs[4*(i-1)+2],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(length(subs)<(4*(i-1)+3)){
      graph3 <- NULL
    } else {
      print(subs[4*(i-1)+3])
      lenGraph(df[df$target==subs[4*(i-1)+3],], seqLen, ntCap, sub=subs[4*(i-1)+3], name=name, target=target)
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+3],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph3 <- plot_grid(title, nrgGraph(df[df$target==subs[4*(i-1)+3],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(length(subs)<(4*(i-1)+4)){
      graph4 <- NULL
    } else {
      print(subs[4*(i-1)+4])
      lenGraph(df[df$target==subs[4*(i-1)+4],], seqLen, ntCap, sub=subs[4*(i-1)+4], name=name, target=target)
      titleText <- gsub("\\.","",paste(subs[4*(i-1)+4],target, sep=" "))
      title <- ggdraw() + draw_label(titleText, fontface='bold')
      graph4 <- plot_grid(title, nrgGraph(df[df$target==subs[4*(i-1)+4],], ntCap, name=paste(name, titleText)),ncol=1, nrow=2, rel_heights = c(0.5,5))
    }
    
    if(i==1){
      graph <- plot_grid(graph1, graph2, graph3, graph4, ncol=4, nrow=1, rel_widths = c(5,5,5,5), labels = c(LETTERS[(4*(i-1)+1):(4*(i-1)+4)], NULL))
    } else {
      graph <- plot_grid(graph, plot_grid(graph1, graph2, graph3, graph4, ncol=4, nrow=1, rel_widths = c(5,5,5,5), labels = c(LETTERS[(4*(i-1)+1):(4*(i-1)+4)], NULL)), rel_widths=c(1,1), rel_heights=c(i-1, 1), ncol=1, nrow=2)
    }
    
    rm(graph1, graph2, graph3, graph4)
    gc()
    
    if(i==ceiling(length(subs)/4)){
      print('plotting')
      save_plot(paste(gsub("\\.","",paste(name, "_by_", target, sep="")),'.png',sep=""), graph, base_height = 5.5*i, base_width = 20, dpi=900)
      rm(graph)
    }
  }
}

setwd("./Helical Stability/3nt")
setwd("./Puerto Rico")
pr8 <- dfPreProcess(pr8, subsa=c(1:5), subsb=c(6:8), converter = nrgConverter, strain="Puerto Rico")
gc()

setwd('..')
setwd("./Hong Kong")
hk <- dfPreProcess(hk, subsa=c(1:6), subsb=c(7:8), converter = nrgConverter, strain="Hong Kong")
gc()

setwd('..')
setwd("./WSN")
wsn <- dfPreProcess(wsn, subsa=c(1:5), subsb=c(6:8), converter = nrgConverter, strain="WSN")
#getNrgGraphs(wsn, 3, target = 'Subunit', seqLen = 9:16)
gc()

setwd('..')
setwd("./Brisbane")
bri <- dfPreProcess(bri, subsa=c(1:5), subsb=c(6:8), converter = nrgConverter, strain="Brisbane")
gc()

setwd("..")
global <- rbind(pr8, hk, wsn, bri)
rm(pr8, hk, wsn, bri)
gc()

#turn off intermediates and final transcription
#global <- global[global$Intermediate==F,] #this option does not change the data meanigfully 

getNrgGraphs(global, 3, target = 'Strain', seqLen = 9:16)

#get statistics
global <- global[global$Realign==0 | global$Realign>=3,]
gc()
strains <- c("Puerto Rico", 'Hong Kong', 'WSN', 'Brisbane')
series <- c('Transcribed', 'Will Realign', 'Realigns')

global$series <- 'Transcribed'
global$series[global$Realign == 3] <- 'Realigns'
global$series[global$Realign > 3] <- 'Will Realign'

means <- vector(length = 12, mode = 'numeric') #means
sds <- vector(length = 12, mode = 'numeric') #standard deviation
p1 <- vector(length = 12, mode = 'numeric') #p versus transcription
p2 <- vector(length = 12, mode = 'numeric') #p versus will realign
p3 <- vector(length = 12, mode = 'numeric') #p versus realigns

for(i in 1:4){
  for(j in 1:3){
    means[(i-1)*3+j] <- mean(global$nt_addition_3[global$Strain==strains[i] & global$series==series[j]])
    sds[(i-1)*3+j] <- sd(global$nt_addition_3[global$Strain==strains[i] & global$series==series[j]])
    p1[(i-1)*3+j] <- as.numeric(t.test(global$nt_addition_3[global$Strain==strains[i] & global$series==series[1]], global$nt_addition_3[global$Strain==strains[i] & global$series==series[j]])$p.value)
    p2[(i-1)*3+j] <- as.numeric(t.test(global$nt_addition_3[global$Strain==strains[i] & global$series==series[2]], global$nt_addition_3[global$Strain==strains[i] & global$series==series[j]])$p.value)
    p3[(i-1)*3+j] <- as.numeric(t.test(global$nt_addition_3[global$Strain==strains[i] & global$series==series[3]], global$nt_addition_3[global$Strain==strains[i] & global$series==series[j]])$p.value)
  }
}

statsTable <- data.frame(Strain = c(rep(strains[1], 3), rep(strains[2], 3), rep(strains[3], 3), rep(strains[4], 3)),
                         Series = rep(series, 4),
                         Mean = means,
                         Standard_Deviation = sds,
                         Lower_CI = 0,
                         Upper_CI = 0,
                         pvalue_vs_Transcribed = p1,
                         pvalue_vs_Will_Realign = p2,
                         pvalue_vs_Realigns = p3,
                         stringsAsFactors = F)

statsTable$Upper_CI <- statsTable$Mean-2*statsTable$Standard_Deviation
statsTable$Lower_CI <- statsTable$Mean+2*statsTable$Standard_Deviation
colnames(statsTable) <- gsub('_', ' ', colnames(statsTable))

write.csv(statsTable, 'Energy Stats.csv', row.names = F)

rm(global)
gc()
