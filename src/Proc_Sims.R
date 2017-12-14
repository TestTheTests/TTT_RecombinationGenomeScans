### KE Lotterhos
### 20171125
### This script is use to process the raw output 
###### from the LOW RECOMBINATION simulations. It produces
###### a table (for each replication simulation)
###### with test statistics and P-values (in columns)
###### for each locus (in rows). 
### The column name format is <method>_<version>_<condition>_<stat>
### The script does one simulation at a time, 
###### so it can be pseudo-parallelized on a cluster.
### Filtering: loci with MAF < 0.01 are removed
### Script assumes libraries are already installed
### This file is in the "/src" folder
### Output written to "/results_proc" folder
### Then concatenated into one large file

### load libraries ####
  library(OutFLANK)
  library(vcfR)
  library(adegenet)
  library(pcadapt)
  library(lfmm)
  library(LEA) 
  library(bigsnpr)
  library(bigstatsr)
  #library(RColorBrewer)
  #library(ggplot2)
  #library(fields)
  library(rehh)

### input arguments ####
  args <- commandArgs(TRUE)
  seed <- args[1]

### read files ####
  folder <- "results/"
    # this script will be called from working directory in this folder
  type <- "RecomLowReg"
  vcf <- read.vcfR(paste( seed, "_", type, "_VCFallsim1.vcf.gz", sep=""))
  ind0 <- read.table(paste( seed,  "_", type,"_outputIndAll.txt", sep=""), header=TRUE)
  muts <- read.table(paste(seed, "_", type,"_outputMuts.txt", sep=""), header=TRUE)
  phen_env <-  read.table(paste( seed, "_", type,"_outputPhenEnv.txt", sep=""), header=TRUE)
  sim <- system(paste("grep recomb ",  seed, "_", type,"_outputSim.txt", sep=""), TRUE)
  H <- read.table(paste(seed, "_", type,"_H.txt", sep=""), header=TRUE)

### Mutation table ####
# Count mutations that contribute at least 1% to genetic variance
  muts$pa2 <- round(muts$selCoef^2*muts$freq*(1-muts$freq),3)
  muts$prop=NA
  muts$prop[muts$type=="m2"] <- muts$pa2[muts$type=="m2"]/sum(muts$pa2[muts$type=="m2"])
  which(duplicated(muts$position))
  muts$count <- FALSE
  muts$count[muts$prop>=0.01] <- TRUE
  muts$count[muts$type=="m4" & muts$freq > 0.05] <- TRUE
  muts$simID <- seed
  muts$simIDpos <- paste(muts$simID, muts$position, sep="_")
  write.table(muts, paste0("results_proc/",seed,"_",type, "_muts_proc.txt"))

### Fix "CHROM" in VCF ####
  # Output "CHROM" in VCF from SliM does not show where the 
  # recombination breakpoints are.
  # Also output is unordered
  ### Replace chromsome 1 with actual chromosome positions ####
  ends=c(0,  seq(50000, 500000, by=50000))
    # linkage groups recombination breakpoints 0.5
    # hard-coded from simulations
  dim(vcf@gt)
  vcf@fix[,"CHROM"] <- NA
  POS <- as.numeric(vcf@fix[,"POS"])
  for (i in 1:(length(ends)-1)){
    cond <- POS >= ends[i] &  POS < ends[i+1]
    print(c(ends[i], ends[i+1], sum(cond)))
    vcf@fix[cond,"CHROM"] = i
  }
  #table(vcf@fix[,"CHROM"])
  my_ord <- order(as.numeric(vcf@fix[,"POS"]))
  vcf2 <- vcf
  vcf2 <- vcf[my_ord,]
  vcf <- vcf2
  rm(vcf2)

  
### Assign individuals to populations ####
  toclust <- ind0[,c("x","y")]
  d <- dist(toclust)
  hc <- hclust(d, method="ward.D")
  k <- 39
  group <- cutree(hc, k=k)
  table(group)
 # label populations by their environment
  ind0$group <- group
  # the id should put them back into order
  ind <- ind0
  ind <- ind[,-which(colnames(ind)=="phenotype2")]
  rem_ind <- which(ind$phenotype1>4 | ind$phenotype1 < -4)
  # remove outlier phenotype
  ind <- ind[-rem_ind,]
  head(ind)
  
  plot(ind$envi, ind$phenotype1)
  (phen_env_corr <- cor.test(ind$phenotype1, ind$envi))
  qplot(ind$x, ind$y, colour=ind$phenotype1)
  qplot(ind$x, ind$y, colour=ind$envi)

  
### Convert VCF to 012 format ####
  # Remove fixed loci (or heterozygote in every individual)
  # Remove loci with He <0.01
  geno <- vcf@gt[,-1] 
  # Character matrix containing the genotypes
  # individuals in columns
  # Remove 1st column, which is
  position <- getPOS(vcf) # Positions in bp
  chromosome <- getCHROM(vcf) # Chromosome information
  G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
  G[geno %in% c("0/0", "0|0")] <- 0
  G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
  G[geno %in% c("1/1", "1|1")] <- 2
  ## Remove fixed loci or all heterozygotes ####
  #dim(G)
  #head(G[1:10,1:10])
  G_sub <- G[,-rem_ind]
  a_freq <- rowSums(G_sub)/(2*ncol(G_sub))
  # total # derived alleles / total # individuals
  #rem = c(which(rowSums(G)==0), which(rowSums(G-2)==0)) ## fixed loci
  rem = which(a_freq < 0.01 | a_freq > 0.99)
  # remove loci with MAF < 0.01
  length(rem)
  position[rem]
  #G_all = G[-rem,] , G_sub
  training <- list(G = G_sub[-rem,], position = position[-rem], 
                   chromosome = chromosome[-rem])
  vcf_filt <- vcf[-rem,-rem_ind]
  #rm(G, position, chromosome)
  #rm(vcf) 
  
### Obtain a quasi-independent set of loci ####
  options(bigstatsr.typecast.warning = FALSE)
  G_coded <- add_code256(big_copy(t(training$G),
                                         type="raw"),
                                code=bigsnpr:::CODE_012)
    # puts it in the raw format and stores likelihood genotype probability
  newpc<-snp_autoSVD(G=G_coded,
                     infos.chr = as.integer(training$chromosome),
                     infos.pos = training$position)
    # this is doing SNP pruning - removing correlated SNPs
    # take snps with highest MAF and correlate snps around it
    # Snps with R^2 > 0.2 are removed
    # the subset is the indexes of the remaining SNPs
    # common error: TridiagEigen: failed to compute all the eigenvalues
      # missing values?
      # sum(is.na(training$G))
      # columns with no variation (zero or very low std)
      # hist(apply(training$G, 1, sd), xlim=c(0,1))
      # I have noticed that this error arises when sd < 0.1
      
  str(newpc)
  plot(newpc)
  which_pruned =  attr(newpc, which="subset")
    # which_pruned are the indexes of the quasi-independent 
    # set of loci that are kept after pruning for LD
  training$G_coded <- G_coded
  training$G_pruned <- training$G[which_pruned,]#
  training$which_pruned <- which_pruned
  rm(which_pruned)
  #str(training)
  
### Final dataframe ####
  final_df <- data.frame(pos=training$position, chrom=training$chromosome,
                         simID=simID)
  dim(final_df)
  final_df$quasi_indep <- FALSE
  final_df$quasi_indep[training$which_pruned] <- TRUE
  
### H12 - Hscan ####
  H$pos <- round(H$x/2)+1
  which_matchH <- match(final_df$pos, H$pos)
  notmatched <- which(is.na(which_matchH))
  #check <- cbind(H$pos[which_matchH], final_df$pos)
  #check[notmatched,]
  final_df$Hscan_v1.3_H12 <- H$H[which_matchH]
  final_df$Hscan_v1.3_H12[final_df$pos%%50000 < 1000 | final_df$pos%%50000 > (50000-1000)] <- NA
    # ignore haplotypes at the end of chromosomes
  plot(final_df$pos, final_df$Hscan_v1.3_H12)
  #plot(final_df$pos, final_df$Hscan_v1.3_H12, xlim=c(90000, 110000))
  #plot_layers(y_head=6000, y_arrows = c(6000, 5000))
  
### PCA loadings and scores if all data is used ####
  pca_all <-pcadapt(training$G,K=3)
  final_df$pca_ALL_PC1_loadings <- pca_all$loadings[,1]
  final_df$pca_ALL_PC2_loadings <- pca_all$loadings[,2]
  final_df$pca_ALL_PC3_loadings <- pca_all$loadings[,3]
  ind$pca_ALL_PC1_scores <- pca_all$scores[,1]
  ind$pca_ALL_PC2_scores <- pca_all$scores[,2]
  
  plot(final_df$pos, final_df$pca_ALL_PC1_loadings)
  plot(final_df$pos, final_df$pca_ALL_PC2_loadings)    
  plot(final_df$pos, final_df$pca_ALL_PC3_loadings) 
  
  qplot(ind$pca_ALL_PC1_scores, ind$pca_ALL_PC2_scores, colour=ind$envi, 
        main="Individual scores without LD pruning") #+ 
  qplot(ind$pca_ALL_PC1_scores, ind$pca_ALL_PC2_scores, colour=ind$phenotype1, 
        main="Individual scores without LD pruning")

    
### PCA loadings if pruned data is used ####
  pca_pruned <-pcadapt(training$G_pruned,K=3)
  final_df$pca_PRUNED_PC1_loadings[training$which_pruned] <- pca_pruned$loadings[,1]
  final_df$pca_PRUNED_PC2_loadings[training$which_pruned] <- pca_pruned$loadings[,2]  
  ind$pca_PRUNED_PC1_scores <-  pca_pruned$scores[,1] 
  ind$pca_PRUNED_PC2_scores <-  pca_pruned$scores[,2] 
  
  plot(final_df$pos, final_df$pca_PRUNED_PC1_loadings)
  plot(final_df$pos, final_df$pca_PRUNED_PC2_loadings)
  
  qplot(ind$pca_PRUNED_PC1_scores, ind$pca_PRUNED_PC2_scores, colour=ind$envi, 
        main="Individual scores WITH LD pruning") #+ 
  qplot(ind$pca_PRUNED_PC1_scores, ind$pca_PRUNED_PC2_scores, colour=ind$phenotype1, 
        main="Individual scores WITH LD pruning")
  
  plot(ind$pca_ALL_PC1_scores, ind$envi)
    abline(lm(ind$envi~ind$pca_ALL_PC1_scores))
  plot(ind$pca_PRUNED_PC1_scores, ind$envi)
    abline(lm(ind$envi~ind$pca_PRUNED_PC1_scores))
  
### PCADAPT all data ####
    final_df$pcadapt_3.1.0_ALL_chisq <- as.numeric(pca_all$chi2.stat)
    final_df$pcadapt_3.1.0_ALL_log10p <- -log10(pca_all$pvalues)
    plot(final_df$pos, final_df$pcadapt_3.1.0_ALL_log10p)
### PCADAPT pruned data ####    
    test <- snp_gc(snp_pcadapt(training$G_coded,U.row = newpc$u[,1]))
    final_df$pcadapt_3.1.0_PRUNED_log10p <- -predict(test,log10=T)
    plot(final_df$pos,final_df$pcadapt_3.1.0_PRUNED_log10p)

### OutFLANK all data ####
    FstDataFrame <- MakeDiploidFSTMat(t(training$G),training$position,
                                      ind$group)
    out_ini <- OutFLANK(FstDataFrame, NumberOfSamples=k) 
    str(out_ini)
    OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
                           NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                             FALSE, RightZoomFraction = 0.05, titletext = NULL)
    OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
                           NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                             TRUE, RightZoomFraction = 0.15, titletext = NULL)
    identical(final_df$pos, FstDataFrame$LocusName)
    final_df$OutFLANK_0.2_FST <- FstDataFrame$FST
    final_df$OutFLANK_0.2_He <- FstDataFrame$He
    
    identical(final_df$pos, final_df$out_ini$results$LocusName)
    head(out_ini$results)
    head(final_df)
    reorder <- match(out_ini$results$LocusName, final_df$pos)
    identical(out_ini$results$LocusName[reorder], final_df$pos)
    final_df$OutFLANK_0.2_ALL_log10p <- -log10(out_ini$results$pvaluesRightTail[reorder])
    hist(out_ini$results$pvaluesRightTail)
    plot(final_df$pos, final_df$OutFLANK_0.2_ALL_log10p)

### OutFLANK pruned data ####    
    out_pruned <- OutFLANK(FstDataFrame[training$which_pruned,], NumberOfSamples=k)     
    str(out_pruned)
    
    P1 <- pOutlierFinderChiSqNoCorr(FstDataFrame, 
                                    Fstbar = out_pruned$FSTNoCorrbar, 
                                    dfInferred = out_pruned$dfInferred, Hmin=0.1)
    P1 <- P1[order(P1$LocusName),]
    identical(P1$LocusName, final_df$pos)
    final_df$OutFLANK_0.2_PRUNED_log10p <- -log10(P1$pvaluesRightTail)
    #plot(final_df$pos, final_df$OutFLANK_0.2_He)
    plot(final_df$pos, final_df$OutFLANK_0.2_PRUNED_log10p) 
    
    plot(final_df$OutFLANK_0.2_ALL_log10p, final_df$OutFLANK_0.2_PRUNED_log10p)

### LFMM Ridge regression model genotype ~ phenotype ####
    # extract scaled genotypes
    scaled.genotype <- scale(as.matrix(t(training$G)))
    #scaled.genotype <- as.matrix(t(sim1$G))
    # extract scaled phenotypes
    phen <- scale(as.matrix(ind$phenotype1))
    # ridge regression
    lfmm.ridge <- lfmm::lfmm_ridge(Y = scaled.genotype, X = phen, K = 3, lambda = 1e-4)
    #The lfmm.ridge object contains estimates for the latent variables and for the effect sizes. Here, the estimates are used for computing calibrated significance values and for testing associations between the response matrix Y and the explanatory variable x. It can be done as follows:
    lfmm.test.ridge <- lfmm::lfmm_test(Y = scaled.genotype, X = phen, lfmm = lfmm.ridge, calibrate = "gif")
    p.values.ridge <- lfmm.test$calibrated.pvalue
    lfmm.test.ridge$gif
    hist(p.values.ridge, col = "lightgreen", main="LFMM ridge")
    #qval <- qvalue::qvalue(p.values)
    #plot(qval)
    #The plot suggests that setting fdr.level = 0.025 warrant few false positives.
    #qval <- qvalue::qvalue(p.values, fdr.level = 0.005)
    #candidates <- which(qval$significant)
    plot(training$position, -log10(p.values.ridge), cex = .5, pch = 19, col = "black", main="LFMM ridge", ylim=c(0, 60))
    plot_layers(y_head=55, y_arrows=c(10, 0))
    
    final_df$LFMM_ridge_0.0_ALL_log10p <- as.numeric(p.values.ridge)
 
### LFMM LASSO model model genotype ~ phenotype ####    
    #LFMM parameters can alternatively be estimated by solving regularized least-squares mimimization, with lasso penalty as follows.
    lfmm.lasso <- lfmm::lfmm_lasso(Y = scaled.genotype, X = phen, K = 3, nozero.prop = 0.02)
    # The lfmm.lasso object contains new estimates for the latent variables and for the effect sizes. 
    # Note that for lasso, we didn't set the value of a regularization parameter. 
    # Instead, we set the proportion of non-null effects (here 2 percent).
    lfmm.test.lasso <- lfmm::lfmm_test(Y = scaled.genotype, X = phen, lfmm = lfmm.lasso, calibrate = "gif")
    p.values.lasso <- lfmm.test.lasso$calibrated.pvalue
    lfmm.test.lasso$gif
    hist(p.values.lasso, col = "lightblue")
    plot(training$position, -log10(p.values.lasso), cex = .5, pch = 19, col = "black", main="LFMM lasso", xaxs="i")
    #plot_layers(y_head=45, y_arrows=c(5,0))
    
    final_df$LFMM_lasso_0.0_ALL_log10p <- as.numeric(p.values.lasso)
    

### Bayesian (LEA) model genotype ~ environment ####
    # Creation of the genotype file: "genotypes.lfmm"
    # 400 SNPs for 50 individuals.
    pc <- prcomp(scaled.genotype)
    plot(pc$sdev[1:10]^2, pch = 19)
    gename <- paste0(seed, "_genotypes.lfmm")
    envname <- paste0( seed, "_gradients.env")
    write.lfmm(t(training$G), gename)
    write.env(ind$envi,envname)
    ################
    # running lfmm #
    ################
    project = lfmm(gename, 
                    envname, 
                    K = 1:3, 
                    repetitions = 3, 
                    project = "new")
    #project = load.lfmmProject(paste0(sub(".lfmm","",gename),"_",seed,"_gradients.lfmmProject"))
    # summary of the project
    summary(project)
    # get adjusted p-values using all runs
    pv1 <-  adjusted.pvalues(project, K = 1)
    pv2 <-  adjusted.pvalues(project, K = 2)
    pv3 <-  adjusted.pvalues(project, K = 3)
    # get the z-scores
    z_ave1 <- rowMeans(z.scores(project, K = 1))
    z_ave2 <- rowMeans(z.scores(project, K = 2))
    z_ave3 <- rowMeans(z.scores(project, K = 3))

    plot(training$position, -log10(pv1$p.values), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA",  xaxs="i")
    plot(training$position, -log10(pv2$p.values), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
    plot(training$position, -log10(pv3$p.values), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
    
    plot(training$position, abs(z_ave1), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
    plot(training$position, abs(z_ave2), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
    plot(training$position, abs(z_ave3), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
    
    final_df$LEA_1.8.1_ALL_K1_log10p <-  -log10(pv1$p.values)
    final_df$LEA_1.8.1_ALL_K2_log10p <-  -log10(pv2$p.values)
    final_df$LEA_1.8.1_ALL_K3_log10p <-  -log10(pv3$p.values)
    final_df$LEA_1.8.1_ALL_K1_z <- z_ave1
    final_df$LEA_1.8.1_ALL_K2_z <- z_ave2
    final_df$LEA_1.8.1_ALL_K3_z <- z_ave3

### Raw associations with environment ####
    spearmans <- apply(training$G, 1, FUN=function(x){
      a <- cor.test(x, ind$envi, method = "spearman")
      return(a$estimate)
    })
    length(spearmans)
    sum(is.na(spearmans))
    hist(spearmans)
    plot(training$position, abs(spearmans), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="raw association", xaxs="i")
    
    final_df$Spearmans_ALL_rho <- spearmans
    
### iHS convert ####
    nlociperchrom <- table(vcf@fix[,"CHROM"])
    
    substring("0|1", 1, last=1)
    substring("0|1", 3, last=3)
    
    head(vcf_filt)
    rem2 <- which(duplicated(vcf_filt@fix[,"POS"]))
    rem2
    
    vcf_filt2 = vcf_filt[-rem2,]
    
    ### Get into right format
    for (i in 1:length(nlociperchrom)){
      keep <- which((vcf_filt2@fix[,"CHROM"]==i))
      #head(vcf_filt2@gt[keep,1:10], 10)
      hap1 <- apply(vcf_filt2@gt[keep,-1], 1, FUN=function(x) substring(x,1,1))
      #dim(hap1)
      #head(hap1[,1:10])
      hap2 <-  apply(vcf_filt2@gt[keep,-1], 1, FUN=function(x) substring(x,3,3))
      #head(hap2[,1:10])
      
      nind <- nrow(ind)
      hapt_out <- matrix(NA, nrow=2*nind, ncol=length(keep))
      odd <- seq(1,(2*nind), by=2)
      even <- odd +1
      hapt_out[odd,] <- as.numeric(hap1)
      rownames(hapt_out) <- rep("", nrow(hapt_out))
      rownames(hapt_out)[odd] <- rownames(hap1)
      rownames(hapt_out)[even] <- rownames(hap2)
      hapt_out[even,] <- as.numeric(hap2)
      #head(hapt_out[,1:10])
      write.table(cbind(rownames(hapt_out), hapt_out+1), paste0("results/",seed,"_chrom",i,".thap"), row.names=F, col.names=F, quote = FALSE)
    }
    #a<- read.table("chrom1.thap")
    
    ### Also need to convert map.inp
    #Each line contains five columns corresponding to:
    #1. the SNP name
    #2. the SNP chromosome (or scaffold) of origin
    #3. the SNP position on the chromosome (or scaffold). Note that it is up to the user to choose either
    #physical or genetic map positions (if available).
    #4. the SNP ancestral allele (as coded in the haplotype input file)
    #5. the SNP derived alleles (as coded in the haplotype input file)
    map <- data.frame(name=1:nrow(vcf_filt2), chrom=as.numeric(vcf_filt2@fix[,"CHROM"]), 
                      pos=as.numeric(vcf_filt2@fix[,"POS"]), anc=1, derived=2)
    # setting anc=0 and derived = 1 thinks missing data
    head(map)
    which(duplicated(map$pos))
    write.table(map, paste0("results/",seed,"_map.inp"), row.names=F, col.names=F)

### iHS run ####
    cnt=0
    wgscan <- NULL
    for(i in 1:length(nlociperchrom)){
      cnt=cnt+1
      tmp.hapfile=paste0("results/",seed,"_chrom",i,".thap")
      
      tmp.hap=data2haplohh(hap_file=tmp.hapfile, map_file=paste0("results/",seed,"_map.inp"), chr.name=i,haplotype.in.columns=FALSE)
      
      tmp.scan=scan_hh(tmp.hap,threads=4)
      
      if(cnt==1){wgscan=tmp.scan}else{wgscan=rbind(wgscan,tmp.scan)}
    }
    
    head(wgscan)
    tail(wgscan)
    
    ihs=ihh2ihs(wgscan,minmaf=0.05,freqbin=0.05)
    head(ihs$iHS,25)
    colnames(ihs$iHS)[4] <- "rehh_2.0.2_ALL_log10p"
    colnames(ihs$iHS)[3] <- "rehh_2.0.2_ALL_iHS"
    #tail(ihs$iHS)
    ihs$frequency.class
    
    plot(ihs$iHS$POSITION,ihs$iHS$`-log10(p-value)`, col=rgb(0,0,0,0.2), pch=20, main="REHH iHS")
    plot(ihs$iHS$POSITION,ihs$iHS$iHS, col=rgb(0,0,0,0.2), pch=20, main="REHH iHS")
 
    final_df <- merge(final_df, ihs$iHS, by.x = "pos", by.y="POSITION")

    head(final_df)    

### Tables to write to file ####
    write.vcf(vcf_filt, file = paste("../results_final/", seed, "_", type, "_VCFallFILT.vcf.gz", sep=""))
    write.table(ind, file = paste("../results_final/", seed, "_", type, "_indFILT.txt", sep=""), row.names=FALSE)
    write.table(final_df, file = paste("../results_final/", seed, "_", type, "_ScanResults.txt", sep=""), row.names = FALSE)
    # sim metrics to write to file
    #simID
    #simtype
    #neutralFSTpruned
    #neutalFSTall
    #neutraldfpruned
    #neutraldfall
    
    