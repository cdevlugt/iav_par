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

templates <- c('CGUU', 'CGUU', 'CGUUU', 'CGUUU', 'CGUUUU', 'CGUUUU', 'CGUUUUC', 'CGUUUUC') #3' vRNA
additions <- c('GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA') #sequence transcribed

templates <- c(templates, 'CGCU', 'CGCU', 'CGCUU', 'CGCUU', 'CGCUUU', 'CGCUUU', 'CGCUUUC', 'CGCUUUC') #3' vRNA
additions <- c(additions, 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA') #sequence transcribed

dfArgs <- list(Additions=additions, Templates=templates, stringsAsFactors=F)

nrgConverter <- do.call(data.frame, dfArgs)
rm(dfArgs)

nrgConverter$TotalEnergy <- helixEnergy(x53=nrgConverter[,1], y35=nrgConverter[,2], offset=0, calc='max')

#setwd() #if not using default
setwd("./Helical Stability/") #folder must be created

#pull GCA and GCG
gca <- nrgConverter[seq(1, nrow(nrgConverter), by = 2),]
gcg <- nrgConverter[seq(2, nrow(nrgConverter), by = 2),]
nrgConverter <- data.frame(template = unique(nrgConverter$Templates),
                           addition = c('GC[A/G]', 'GC[A/G]A', 'GC[A/G]AA', 'GC[A/G]AAA'),
                           GCA = gca$TotalEnergy,
                           GCG = gcg$TotalEnergy,
                           stringsAsFactors = F)
nrgConverter$difference <- abs(nrgConverter$GCA - nrgConverter$GCG)
write.csv(nrgConverter, 'Energy Additions.csv', row.names = F)