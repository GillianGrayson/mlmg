################################################################################
# Canadian Bioinformatics Workshop 2017
# Bioinformatics for Genomic Medicine: Module 6
#
# Created:  11-May-2017.
# Author: Andrei Turinsky
###############################################################################

source("http://bioconductor.org/biocLite.R");

biocLite(FlowSorted.Blood.450k);
library(GEOquery);
library(gplots);
library(lumi);
biocLite(minfi);
library(matrixStats);
biocLite(minfiData);
biocLite(sva);

dataPath = '.'

##########################################################
# Load data matrix

gse_list = getGEO(GEO = "GSE52588", GSEMatrix=TRUE, getGPL=FALSE)

class(gse_list)
length(gse_list)

gse_eset = gse_list[[1]]
class(gse_eset)

# Alternatively:
# gse_eset = getGEO(
#		filename= file.path(dataPath, 'GSE52588_series_matrix.txt.gz'),
#		GSEMatrix=TRUE, getGPL=FALSE)


##########################################################
# Load Illumina HumanMethylation450 array annotation: pre-load from
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL13534&form=text&view=full

gpl = getGEO( filename=file.path(dataPath, 'GPL13534.txt'))

annotationTable = Table(gpl)
head(annotationTable)
rownames(annotationTable) = annotationTable$ID

##########################################################
# Extract phenotype data	

pheno = pData(gse_eset)

head(pheno)
pheno = pheno[, grep("characteristics_ch1", colnames(pheno))]
head(pheno)

levels(pheno$characteristics_ch1)

groups = as.character(pheno$characteristics_ch1)
groups = gsub("su", "ds", substr(groups, 8,9))

family = as.character(pheno$characteristics_ch1.2)
family = gsub("family: ", "", family)

pheno$Sample_Name  = paste(groups, family, sep="")
pheno$Sample_Group  = groups
pheno$gsm = rownames(pheno)

head(pheno)

##########################################################
# Create a minfi methylation data object	

beta_matrix = exprs(gse_eset)
head(beta_matrix)

# exclude "mother" samples

isNotMother = pheno$Sample_Group != "mo"
isNotMother

beta_matrix = beta_matrix[, isNotMother]
pheno = pheno[isNotMother,]

# remove probes with missing and extreme values

hasValueNA = rowSums( is.na(beta_matrix) ) > 0
hasValue0  = rowSums( beta_matrix == 0 , na.rm = TRUE) > 0  
hasValue1  = rowSums( beta_matrix == 1 , na.rm = TRUE) > 0 

# remove probes on sex chromosomes

probes = rownames(beta_matrix)
isChrXY = as.character(annotationTable[ probes, 'CHR']) %in% c('X', 'Y')
isNonCG = substr(probes, 1, 2) != 'cg'

cat( sum(hasValueNA), "probes have too many missing values\n")
cat( sum(hasValue0),  "probes have too many values = 0\n")
cat( sum(hasValue1),  "probes have too many values = 1\n")
cat( sum(isChrXY),    "probes are on sex chromosomes\n")
cat( sum(isNonCG),    "probes are on not CpG\n")

isGoodProbe = !hasValueNA & !hasValue0 & !hasValue1 & !isChrXY & !isNonCG
beta_matrix = beta_matrix[isGoodProbe,]

dim(beta_matrix)

colnames(beta_matrix) = pheno$Sample_Name
rownames(pheno) = pheno$Sample_Name

# Make a GenomeRatioSet object, see https://support.bioconductor.org/p/73941

require(minfiData)

tmpRSet = RatioSet(Beta = beta_matrix) 
annotation(tmpRSet) = annotation(MsetEx) 
pData(tmpRSet) = pheno

grSet = mapToGenome(tmpRSet)	

# Note: can also do: annotationTable = getAnnotation(MsetEx) 
# but beware of slightly different column names, some missing probes

##########################################################
# Find disease signature

beta_matrix = getBeta(grSet)
m_matrix = getM(grSet)
pheno = pData(grSet)

require(lumi)
m2beta
beta2m

samplesDisease = rownames(pheno)[ pheno$Sample_Group == 'ds']
samplesControl = rownames(pheno)[ pheno$Sample_Group == 'si']

# Find differentially methylated positions

dmp = dmpFinder(m_matrix, pheno$Sample_Group, type = "categorical")
head(dmp)

adjustMethod = 'bonf'

dmpCpgs = rownames(dmp)
dmp$betaDisease = rowMeans( beta_matrix[dmpCpgs, samplesDisease, drop=F]) 
dmp$betaControl = rowMeans( beta_matrix[dmpCpgs, samplesControl, drop=F]) 
dmp$deltaBeta = dmp$betaDisease - dmp$betaControl
dmp$adjPval = p.adjust(dmp$pval, method = adjustMethod)
dmp$genes = annotationTable[dmpCpgs, "UCSC_RefGene_Name"]
dmp$chr = annotationTable[dmpCpgs, "CHR"]

head(dmp)

# extract signature

pvalThresh = 0.01
deltaBetaThresh = 0.20

sigCpgs = dmpCpgs[ (dmp$adjPval < pvalThresh) & (abs(dmp$deltaBeta) > deltaBetaThresh) ]

cat("Found", length(sigCpgs), "signature CpGs:", "\n")
cat("\t", sum((dmp$pval < pvalThresh), na.rm=T), "CpGs with pval <", pvalThresh, "\n")
cat("\t", sum((dmp$adjPval < pvalThresh), na.rm=T), "CpGs with pval <", pvalThresh, "after", adjustMethod, "adjustment\n")
cat("\t", sum(abs(dmp$deltaBeta) > deltaBetaThresh, na.rm=T), "CpGs with abs(deltaBeta) >", deltaBetaThresh, "\n")
cat("\t", sum((dmp$adjPval < pvalThresh) & (abs(dmp$deltaBeta) > deltaBetaThresh), na.rm=T), "CpGs satisfy both\n")

# check signature CpG chromosomes

summary(annotationTable[sigCpgs,'CHR'])

##########################################################
# hierarchical clustering

beta_matrix = getBeta(grSet)
beta_matrix_sig = beta_matrix[sigCpgs,]


dist_global = dist( t(beta_matrix) )
hc_global = hclust(dist_global)
plot(hc_global)

dist_signature = dist( t(beta_matrix_sig))
hc_signature = hclust(dist_signature)
plot(hc_signature)

plot( as.dendrogram(hc_signature))

hc_signature = hclust(dist_signature, method="ward.D2")
plot(hc_signature)

dist_signature = dist( t(beta_matrix_sig), method = 'manhattan')
hc_signature = hclust(dist_signature, method="ward.D2")
plot(hc_signature)

# extract cluster membership

clusters_signature = cutree(hc_signature, k = 2)	
clusters_signature

##########################################################
# heat maps

# basic heatmap

heatmap(beta_matrix_sig)

# better heatmaps

require(gplots)

heatmap.2(beta_matrix_sig)
heatmap.2(beta_matrix_sig, trace="none")

# make a pretty heatmap

dendrogram = "column" # c("both","row","column","none")
method = 'manhattan' 
clustMethod = "ward.D2" 

colorRampPalette = colorRampPalette(c("blue", "white",  "orange"))(n = 100)
sampleColors = ifelse(pheno$Sample_Group == 'ds', 'red', 'turquoise')

symmetric = FALSE
reorder = TRUE


heatmap.2( 
  beta_matrix_sig, 
  dendrogram = dendrogram,
  Rowv = reorder,
  Colv=if(symmetric) "Rowv" else reorder,
  symm=symmetric, 
  symkey=F, 
  symbreaks=F,
  trace="none", 
  density.info="density", #"density" "none"
  distfun = function(x) { dist(x, method = method) },
  hclustfun = function(x) hclust(x, method = clustMethod),
  col =  colorRampPalette, 
  key.title = "DNAm scale", 
  ColSideColors = sampleColors,
  labRow = ''
)

##########################################################
# Principal component analysis plots

pca_global = prcomp(t(beta_matrix))
pca_signature = prcomp(t(beta_matrix_sig))

samples = rownames(pheno)
sampleColors = ifelse(pheno$Sample_Group == 'ds', 'red', 'turquoise')

par(mfrow = c(1,2))

plot(pca_global$x[,1:2], type='n', main="PCA - global")
text(pca_global$x[,1:2], samples, col = sampleColors)

plot(pca_signature$x[,1:2], type='n', main="PCA - signature")
text(pca_signature$x[,1:2], samples, col = sampleColors)

par(mfrow = c(1,1))


##########################################################
# Influence of cell type composition

require(FlowSorted.Blood.450k)
data(FlowSorted.Blood.450k.compTable)
data(FlowSorted.Blood.450k.JaffeModelPars)

head(FlowSorted.Blood.450k.JaffeModelPars)

isSuspicious = p.adjust(FlowSorted.Blood.450k.compTable$p.value, 'bonf') < 0.05
suspiciousBloodCpgs_all = rownames(FlowSorted.Blood.450k.compTable)[ isSuspicious ]

mean( rownames(grSet) %in% suspiciousBloodCpgs_all)
mean( sigCpgs %in% suspiciousBloodCpgs_all)

sigCpgs_clean = setdiff(sigCpgs, suspiciousBloodCpgs_all)
sigCpgs_clean

# Top 600 CpGs variable used in Jaffe-Irizarry cell subtype model

suspiciousBloodCpgs_top600 = rownames(FlowSorted.Blood.450k.JaffeModelPars)

mean(rownames(grSet) %in% suspiciousBloodCpgs_top600)
mean(sigCpgs %in% suspiciousBloodCpgs_top600)


##########################################################
# Classification plot - similar to (Choufani et al., 2015)

# split data into training and test sets, find a training-set signature

set.seed(1234)
trainSamples = sample( colnames(grSet), size = 20)
testSamples = setdiff(colnames(grSet), trainSamples)

dmp_train = dmpFinder(m_matrix[,trainSamples], pheno[trainSamples,]$Sample_Group, type = "categorical")
head(dmp_train)

train_sigCpgs = rownames(dmp_train)[ dmp_train$qval < 0.01]
train_sigCpgs = setdiff(train_sigCpgs, suspiciousBloodCpgs_all)

# build DNAm profiles and score the training set

trainingDisease = intersect(trainSamples, samplesDisease)
trainingControl = intersect(trainSamples, samplesControl)

sigProfileDisease = rowMedians(beta_matrix[train_sigCpgs, trainingDisease])
sigProfileControl = rowMedians(beta_matrix[train_sigCpgs, trainingControl])

coords = sapply(testSamples, function(newSample) {
  sigProfile = beta_matrix[train_sigCpgs, newSample]
  x = cor(sigProfile, sigProfileControl)
  y = cor(sigProfile, sigProfileDisease)
  c(x,y)
} )
coords = t(coords)

# plot the classification score

minCorr = .8
plot( NULL, xlim=c(minCorr,1), ylim=c(minCorr,1), 
      xlab = "Similarity to control", ylab = "Similarity to disease")
abline(a=0, b=1, col = 'pink')

colors = ifelse(pheno[testSamples, "Sample_Group"] == 'ds', 'red', 'turquoise')
text(coords, testSamples, col = colors)


##########################################################
# Batch correction

require(minfiData)

stopifnot( all(rownames(grSet) %in% rownames(MsetEx)))

# simulate a contaminant: use the last sample of MsetEx dataset

NUM_BATCH_CPGS = 10000

batchCpgs = rownames(grSet)[1:NUM_BATCH_CPGS]
contaminant = getBeta(MsetEx)[batchCpgs, 6]

# add 15-20% contamination to randomly chosen samples

set.seed(1234)
compromisedSamples = sample( colnames(grSet), size = 30)

batch_beta = getBeta(grSet)[batchCpgs,] 
for(j in compromisedSamples) {
  q = runif(n=1, min = .15, max = .20)
  batch_beta[, j] = q * contaminant + (1-q) * batch_beta[,j]
}

# PCA plot

sampleNames = colnames(batch_beta)
isCompromized = (sampleNames %in% compromisedSamples)

pca_batch_compr = prcomp(t(batch_beta))
plot(pca_batch_compr$x[,1:2], type='n', main="PCA - batch")
text(pca_batch_compr$x[,1:2], sampleNames, col = ifelse(isCompromized, 'red', 'blue'))

# deviations of contaminated data from the original data

boxplot( getBeta(grSet)[batchCpgs,] - batch_beta, col = 'grey', las=2)

# batch correction

require(sva)

batch = as.numeric(isCompromized)
modcombat = model.matrix(~ Sample_Group, data=pheno)
combat_beta = ComBat(dat=batch_beta, batch=batch, mod=modcombat)


# deviations of batch-corrected data from the original data

boxplot( getBeta(grSet)[batchCpgs,] - combat_beta,  col = 'grey', las=2)

##########################################################
# DMRs Bumps

designMatrix = model.matrix(~ Sample_Group, data = pheno)

# Explore the number of DMRs with different cutoffs, without permutations 

clean_probes = setdiff( rownames(grSet), suspiciousBloodCpgs_all)

for( cutoff in seq(.1, .5, by = .1)) {
  print(cutoff)
  dmrs = bumphunter(grSet[clean_probes,], design = designMatrix, cutoff = cutoff, B=0, type="Beta")
}

# run the bump hunting

date()
dmrs = bumphunter(grSet[clean_probes,], design = designMatrix, cutoff = .3, B=10, type="Beta")
date()

head(dmrs$table)

selected_dmr_table = subset(dmrs$table, L >= 3) 

# plot DMR -  similar to Bioconductor package DMRcate

plotDMR = function(selectedDMR, grSet) {
  
  dataTable = getBeta(grSet)
  
  # retrieve array annotations
  
  arrayCpgs = rownames(annotationTable)
  chromosomes = paste0('chr',annotationTable$CHR)
  positions = annotationTable$MAPINFO
  names(positions) = arrayCpgs
  
  # find matching array probes
  
  chr = selectedDMR[,'chr'] 
  start = selectedDMR[,'start']
  end = selectedDMR[,'end']
  
  isMatchingCpg =  (chromosomes == chr) & (positions >= start) & (positions <= end)
  
  # match to the CpGs in the data and order by chromosomal position
  
  matchingCpgs = arrayCpgs [isMatchingCpg]
  matchingCpgs = intersect(matchingCpgs, rownames(dataTable))
  matchingCpgs = matchingCpgs[ order(positions[matchingCpgs]) ]
  
  matchingPositions = positions[matchingCpgs]
  
  # get a list of matching genes
  
  genes = as.character(annotationTable[matchingCpgs,'UCSC_RefGene_Name'])
  geneName = paste(sort(unique(unlist(strsplit(genes, ';')))), collapse = ', ')
  
  cat("Found ", length(matchingCpgs), "CpGs in the DMR\n")
  cat("REGION:", chr, start, end, '\n', "\tCPGS:", matchingCpgs, '\n')
  cat("\tGenes: ", geneName)
  
  # plot the DMR
  
  colorDisease = rgb(1,0,0, .5)
  colorControl = rgb(0,0,1, .5)
  
  plot(NULL, xlim=c(start, end), ylim=c(0,1), xlab  = paste0("DMR on ", chr, ": ", geneName) , ylab = "DNAm")
  
  # plot each CpG
  
  for(i in 1:length(matchingCpgs)) {
    cpg = matchingCpgs[i]
    pos = matchingPositions[i]
    points( cbind( pos, dataTable[cpg, samplesDisease]), pch=19, col=colorDisease)
    points( cbind( pos, dataTable[cpg, samplesControl]), pch=4,  col=colorControl)
  }
  
  # plot averages as lines
  
  lines(matchingPositions, rowMeans(dataTable[matchingCpgs, samplesDisease]), col=colorDisease)
  lines(matchingPositions, rowMeans(dataTable[matchingCpgs, samplesControl]), col=colorControl)
}

plotDMR( selected_dmr_table[1,], grSet)
plotDMR( selected_dmr_table[2,], grSet)

##########################################################
# End of script