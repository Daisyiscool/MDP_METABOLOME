setwd('/home/daisy/lab_work/natasha')
library(dplyr)
human.full<-read.csv("human_compoundlist_compound_name.csv", sep=",", na.strings=NA, header=F)

human.full =human.full[-1,]
human.full[1:5, 1:5]
a = human.full[1,]
colnames(human.full) = human.full[1,]
human.full =human.full[-1,]
human.full[1:5, 1:5]
#rownames(human.full) = human.full[,1]
human.full <- human.full%>%
  select(ID, everything())
#stool_neg <- stool_neg[, c(1,3:61)]
#stool_neg$Match <- ifelse(is.na(stool_neg$Match), stool_neg$Matched.Compound, stool_neg$Match)
#stool_neg <- stool_neg[, c(1,3:61)]
#stool_neg <- stool_neg %>%
#  distinct(Match, IR_S10, .keep_all = TRUE)
#stool_neg <- stool_neg %>%
#  distinct(Match, .keep_all = TRUE)
#neg <- stool_neg

human.full[,2:17] <- as.data.frame(lapply(human.full[,2:17], as.numeric))#convert all charactor to numeric, R recognize my table as all string not numeric
human<- human.full %>%
  mutate(S1002_01_norm = S1002_01/median(S1002_01, na.rm=TRUE),
         S1004_01_norm = S1004_01/ median(S1004_01, na.rm=TRUE),
         S1006_01_norm = S1006_01 / median(S1006_01, na.rm=TRUE),
         S1013_01_norm = S1013_01 / median(S1013_01, na.rm=TRUE),
         S1014_01_norm = S1014_01 / median(S1014_01, na.rm=TRUE),
         S1020_01_norm = S1020_01 / median(S1020_01, na.rm=TRUE),
         S1028_01_norm = S1028_01 / median(S1028_01, na.rm=TRUE),
         S1037_01_norm = S1037_01 / median(S1037_01, na.rm=TRUE),
         S2006_01_norm = S2006_01 / median(S2006_01, na.rm=TRUE),
         S2007_01_norm = S2007_01 / median(S2007_01, na.rm=TRUE), 
         S2011_01_norm = S2011_01 / median(S2011_01, na.rm=TRUE), 
         S2015_01_norm = S2015_01 / median(S2015_01, na.rm=TRUE),
         S2017_01_norm = S2017_01 / median(S2017_01, na.rm=TRUE),
         S2029_01_norm = S2029_01 / median(S2029_01, na.rm=TRUE), 
         S2036_01_norm = S2036_01 / median(S2036_01, na.rm=TRUE), 
         S2041_01_norm = S2041_01 / median(S2041_01, na.rm=TRUE))
library(tidyverse)
human <- human %>% 
  gather(Sample, Intensity, 18:33)

human <- human[c(1, 18,19)]
human$Sample <-gsub('_norm', '', human$Sample)

#add metadata
human$Group <- ifelse(human$Sample %in% c('S2006_01', 'S2007_01', 'S2011_01', 'S2015_01', 'S2017_01', 'S2029_01','S2036_01', 'S2041_01'), 
                      'Control', 'MDP')
RS =read.csv('human_compoundlist_RS_NR.csv')
human$Responder <- ifelse(human$Sample %in% c('S2006_01', 'S2007_01', 'S2011_01', 'S2015_01', 'S2017_01', 'S2029_01','S2036_01', 'S2041_01'), 'Control', 
                               ifelse(human$Sample %in% c('S1002_01', 'S1004_01', 'S1028_01'), 'Non-responder', 'Responder'))
human_pca<-human%>%
  spread(ID, Intensity)#convert to wideformat
#remove columns with na's 
human_pca<-human_pca%>%
  select_if(~ !any(is.na(.)))

write.csv(human_pca, 'human.compound.normalized_peak_intensity.csv')
MDP_pca <- human_pca %>%
  filter(Group== 'MDP')
write.csv(MDP_pca, 'diet_neg_normalized_peak_intensity.csv')

####now we add log transformed value to the main long table--do not do this before this step
human$log_norm<-log10(1+human$Intensity)
#only keep log_norm column, delete intensity column
log_pca <- human%>% select(-Intensity)
log_pca<-log_pca%>%
  spread(ID, log_norm)
write.csv(log_pca, 'log_transformed_peak_intensity_human.csv')
MDP_logpca <- log_pca %>%
  filter(Group== 'MDP')
write.csv(MDP_logpca, 'MDP_log_transformed_peak_intensity.csv')

###this scaled pca uses intensity value
###human_pca contains intensity value
###log_pca contains log transformed value
scaled.pca<-prcomp(human_pca[c(4:702)],scale=TRUE, center=TRUE) 
scaled.MDP.pca<-prcomp(MDP_pca[c(4:702)],scale=TRUE, center=TRUE) 
library(ggfortify)
responder_cols<-c(
                  "#DFC27D", #non-responder
                  "#80CDC1")
MDP_pca$Responder <- factor(MDP_pca$Responder)
pcaPlot1<-ggplot2::autoplot(scaled.MDP.pca, data=MDP_pca, frame.colour="Responder", loadings=FALSE, colour="Responder", 
                            loadings.label.colour="black", loadings.colour="black", loadings.label=FALSE, frame=TRUE, 
                            loadings.label.size=5, loadings.label.vjust=-1, size=5) + 
  scale_color_manual(values = responder_cols)+
  scale_fill_manual(values = responder_cols)+
  theme_classic()+
  theme(legend.text = element_text(size=18), 
        legend.position="right",
        plot.background = element_blank(),
        legend.title = element_text(size=18, face="bold"), 
        axis.text = element_text(size=18, colour = 'black'), 
        axis.title = element_text(size=18,  face="bold"));pcaPlot1                             
ggsave('pca.all.group.png', dpi = 300, width = 10, height = 8)
#############################################################################
#############################################################################
#############################################################################
#plsda calculation
library("RVAideMemoire")
library(mixOmics)
#assigning datasets 
write.csv(human, 'input for X.csv')
human <- read.csv('input for X.csv', header = T)
X <- human%>%
  dplyr::select(Sample, Group, Responder, ID, log_norm)%>%
  spread(ID, log_norm) #generate spread data frame
levels(as.factor(X$Responder))
library(forcats)
X$Responder<-fct_relevel(X$Responder, "Control", "Non-responder", "Responder")
levels(as.factor(X$Responder))
X$Group<-fct_relevel(X$Group, "Control", "MDP")
levels(as.factor(X$Group))
Y <- as.factor(X$Responder) #select treatment names
Y
Y1 = as.factor(X$Group)
Y1
X1<-X[4:702] #pull only data columns
# run PLSDA 
MyResult.plsda <- plsda(X1,Y1) # 1 Run the method
plotIndiv(MyResult.plsda)    # 2 Plot the samples
plotVar(MyResult.plsda, cutoff = 0.6)
levels(Y)
responder_cols<-c("#8C510A", #Control
                  "#DFC27D", #non-responder
                  "#80CDC1") #responder
pch.input = c(1, 5, 19)#select shape
pdf("PLSDA.allgroups.pdf", width=9, height=6)
plotIndiv(MyResult.plsda, col=responder_cols, ind.names = FALSE, legend=TRUE, legend.title = "Groups", 
          ellipse = T, title="plsDA Plot", pch=19, style = "graphics", centroid=FALSE, point.lwd = 2, cex=2)
dev.off() 
####this is to do the plsda only between control and MDP
MyResult.plsda.MDP_CT = plsda(X1, Y1)
plotIndiv(MyResult.plsda.MDP_CT)    # 2 Plot the samples
plotVar(MyResult.plsda.MDP_CT, cutoff = 0.7)
levels(Y1)
#View the metabolites that are most highly differentiating metabolites.
#it's better to do VIP between two groups, so subsetting to do Control vs MDP and RS vs NR
MyResult.plsda_responder <- plsda(X1,Y, ncomp=9) #number of components is classes-1
library(RVAideMemoire)
allgroups_VIP <- PLSDA.VIP(MyResult.plsda_responder)
allgroups_VIP_DF <- as.data.frame(allgroups_VIP[["tab"]])
allgroups_VIP_DF
# Converting row names to column
allgroups_VIP_table <- rownames_to_column(allgroups_VIP_DF, var = "Metabolite")
#filter for VIP > 1
responder_VIP_1 <- allgroups_VIP_table %>% 
  filter(VIP >= 1.25)
head(responder_VIP_1)
#plot
pdf("VIP_cutoff_1.25.pdf", width=8, height=11)
responder_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("VIP Score") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off() 
#####################################################################
#####################################################################
#####################################################################
#this is the comparision between Contro and MDP only
MyResult.plsda_group<- plsda(X1,Y1, ncomp=9) 
two.group_VIP <- PLSDA.VIP(MyResult.plsda_group)
two.group_VIP_DF <- as.data.frame(two.group_VIP[["tab"]])
two.group_VIP_DF
# Converting row names to column
two.group_VIP_table <- rownames_to_column(two.group_VIP_DF, var = "Metabolite")
#filter for VIP > 1
group_VIP_1 <- two.group_VIP_table %>% 
  filter(VIP >= 1.25)
head(group_VIP_1)
write_csv(group_VIP_1, "Two.groups_VIPs.csv")

#plot
pdf("two.groups.VIP_cutoff_1.25.pdf", width=8, height=19)
group_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("VIP Score (MDP vs Control)") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off() 
###########################################################################
###########################################################################
####next, we will do the same VIP plot between responders and nonresponders
X.MDP = X%>%
  filter(Group == 'MDP')
head(X.MDP)
Y <- as.factor(X.MDP$Responder) #select treatment names
X.MDP<-X.MDP[4:702] #pull only data columns
# run PLSDA 
MyResult.plsda <- plsda(X.MDP,Y) # 1 Run the method
plotIndiv(MyResult.plsda)    # 2 Plot the samples
plotVar(MyResult.plsda, cutoff = 0.8)
levels(Y)
MDP.responder_cols<-c("#DFC27D", #non-responder
                  "#80CDC1") #responder
pch.input = c(5, 19)#select shape
pdf("PLSDA.pdf", width=9, height=6)
plotIndiv(MyResult.plsda, col=MDP.responder_cols, ind.names = FALSE, legend=TRUE, legend.title = "Groups", 
          ellipse = T, title="plsDA Plot", pch=19, style = "graphics", centroid=FALSE, point.lwd = 2, cex=2)
dev.off() 
MyResult.plsda_MDP <- plsda(X.MDP,Y, ncomp=8) #number of components is classes-1
library(RVAideMemoire)
responder.groups_VIP <- PLSDA.VIP(MyResult.plsda_MDP)
MyResult.plsda_MDP_DF <- as.data.frame(responder.groups_VIP[["tab"]])
MyResult.plsda_MDP_DF
# Converting row names to column
responder_VIP_table <- rownames_to_column(MyResult.plsda_MDP_DF, var = "Metabolite")
#filter for VIP > 1
responder_VIP_1 <- responder_VIP_table %>% 
  filter(VIP >= 1.25)
head(responder_VIP_1)
#plot
pdf("VIP_responder.nonresponder_1.25.pdf", width=8, height=15)
responder_VIP_1 %>%
  arrange(VIP) %>%
  ggplot( aes(x = VIP, y = reorder(Metabolite,VIP,sum))) +
  geom_point() +
  ylab("Metabolite") +
  xlab("VIP Score") +
  ggtitle("VIP Score (Responders vs Non-responders)") +
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
#WGCNA prep
metadata<-log_pca%>%
  dplyr::select(Sample, Group, Responder)
#metadata <- metadata %>%
#  distinct(Patient, ID, .keep_all = TRUE)#reduce to just one set of metadata 
human$ID <- paste('C', human$ID, sep = '_')
human_WGCNA<-human%>%
  mutate(tag=paste(Sample, Responder, sep="_"))%>%
  dplyr::select(tag, ID, log_norm)
MDP_WGCNA = human%>%filter(Group == 'MDP')
MDP_WGCNA = MDP_WGCNA%>%
  mutate(tag = paste(Sample, Responder, sep = '_'))%>%
  dplyr::select(tag, ID, log_norm)
#convert data to wide format
human_WGCNA<-human_WGCNA%>%
  spread(tag, log_norm)
head(human_WGCNA)
rownames(human_WGCNA)<-human_WGCNA$ID
human_WGCNA<-human_WGCNA%>%
  dplyr::select(-ID)

MDP_WGCNA<-MDP_WGCNA%>%
  spread(tag, log_norm)
head(MDP_WGCNA)
rownames(MDP_WGCNA)<-MDP_WGCNA$ID
MDP_WGCNA<-MDP_WGCNA%>%
  dplyr::select(-ID)
#Check that there are no metabolites with 0 counts for all samples. Should return TRUE.  
rowSums(dplyr::count(human_WGCNA)) > 0
rowSums(dplyr::count(MDP_WGCNA)) > 0
library(WGCNA)
library(genefilter)
filt <- filterfun(pOverA(0.5,0.01))
#create filter for the counts data
gfilt <- genefilter(MDP_WGCNA, filt)
#identify genes to keep by count filter
keep <- MDP_WGCNA[gfilt,]
#identify gene lists
n.keep <- rownames(keep)
#gene count data filtered in PoverA, P percent of the samples have counts over A
data_filt <- as.data.frame(MDP_WGCNA[which(rownames(MDP_WGCNA) %in% n.keep),])
#How many rows do we have before and after filtering?
nrow(MDP_WGCNA) #Before
nrow(MDP_WGCNA) #After
####we did not end up filter out any compound since they all meet the criterion, Filtering does not remove any metabolites, all 182 metabolites are used in analysis.
#Checking that all row and column names match. Should return "TRUE"
all(rownames(metadata$tag) %in% colnames(data_filt))
all(rownames(metadata$tag) == colnames(data_filt)) 
## Outlier detection  
sampleTree = hclust(dist(data_filt), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf("outliers_compound.MDP.pdf", width = 35, height = 35)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 0.5, cex.axis = 0.5, cex.main = 0.5)
dev.off()#L-phenylalanine is the outlier
#Transpose such that samples are in rows and metabolites are in columns.  
data_filt <- t(data_filt)
#Look for outlier samples by examining tree of samples  
sampleTree = hclust(dist(data_filt), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# 2. Network construction and consensus module detection  
## Choosing a soft-thresholding power: Analysis of network topology β  
#The soft thresholding power (β) is the number to which the co-expression similarity is raised to calculate adjacency. The function pickSoftThreshold performs a network topology analysis. The user chooses a set of candidate powers, however the default parameters are suitable values.  
# # Choose a set of soft-thresholding powers
allowWGCNAThreads()
powers <- c(c(1:10), seq(from = 2, to=20, by=2)) #Create a string of numbers from 1 through 10, and even numbers from 10 through 20
# 
# # Call the network topology analysis function
sft <-pickSoftThreshold(data_filt, powerVector = powers, verbose = 5)

#Plot the results.  
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# # # Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# # # this line corresponds to using an R^2 cut-off
abline(h=0.6,col="red")
# # # Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#From this data, it appears that our **soft thresholding power is 12**. 
picked_power = 14
temp_cor <- cor       
cor <- WGCNA::cor                                             # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk10 <- blockwiseModules(data_filt,                         # <= input here
                          
                          # == Adjacency Function ==
                          power = 14,               # <= power here
                          networkType = "unsigned",
                          
                          # == Tree and Block Options ==
                          deepSplit = 0,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 5,  
                          
                          maxBlockSize = 1000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time) but it doesn't save a file
                          saveTOMs = F,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)




cor <- temp_cor     # Return cor function to original namespace

# Identify labels as numbers 
mergedColors = netwk10$colors
# Plot the dendrogram and the module colors underneath
unique(netwk10$colors)
pdf("blockwise_module_colors.MDPgroup.pdf")
plotDendroAndColors(
  netwk10$dendrograms[[1]],
  mergedColors[netwk10$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
dev.off()

table(mergedColors)
## Relate modules to sample information  
module_df <- data.frame(
  Metabolite = names(netwk10$colors),
  colors = netwk10$colors
  #colors = labels2colors(netwk$colors)
)
module_df[1:5,]

write.csv(module_df, "metabolite_modules.MDP.only.csv")

# Get Module Eigengenes per cluster
MEs <- moduleEigengenes(data_filt, mergedColors)$eigengenes
# Calculate dissimilarity of module eigengenes 
MEDiss = 1-cor(MEs); 
# Cluster module eigengenes 
METree = hclust(as.dist(MEDiss), method = "average") 
# Plot the result 
#sizeGrWindow(7, 6) 
pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6) 
plot(METree, main = "Clustering of module eigengenes", 
     xlab = "", sub = "") 
MEDissThres = 0.6######剪切高度可修改 
# Plot the cut line into the dendrogram 
abline(h=MEDissThres, col = "red") 
dev.off()

#重新绘制合并后的模块层次聚类图

# Call an automatic merging function 
merge = mergeCloseModules(data_filt, mergedColors, cutHeight = MEDissThres, verbose = 3) 
# The merged module colors 
new.Colors = merge$colors
# Eigengenes of the new merged modules: 
mergedMEs = merge$newMEs 
table(new.Colors) 

#sizeGrWindow(12, 9) 
pdf(file="7_merged dynamic.pdf", width = 9, height = 6) 
plotDendroAndColors(netwk10$dendrograms[[1]], cbind(mergedColors, merged.Colors), 
                    c("Block_1", "Block_2"), 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05) 
dev.off()

# Reorder modules so similar modules are next to each other
MEs <- orderMEs(MEs)
module_order = names(MEs) %>% gsub("ME","", .)
mergedMEs = orderMEs(mergedMEs)
new.module_order = names(mergedMEs) %>% gsub("ME","", .)
# Add Sample names
MEs0 <- MEs
MEs0$ID = row.names(MEs)
new.MEs0 = mergedMEs
new.MEs0$ID = row.names(mergedMEs)
# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-ID) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
new.mME = new.MEs0 %>%
  pivot_longer(-ID) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
mME$ID <- factor(mME$ID, levels = c('S1002_01_Non-responder', 'S1004_01_Non-responder', 'S1028_01_Non-responder', 'S1006_01_Responder',
                                            'S1013_01_Responder', 'S1014_01_Responder', 'S1020_01_Responder', 'S1037_01_Responder'))
ggplot(mME, aes(x=ID, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-Sample Relationships", y = "Modules", fill="corr")
ggsave('block1_module vs samples.png', dpi = 300, width = 8, height = 8)
ggplot(new.mME, aes(x=ID, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18, face = 'bold')) +
  labs(title = "Module-Sample Relationships", y = "Modules", fill="corr")
ggsave('sample.png', dpi = 300, width = 8, height = 8)
unique(new.mME$ID)
new.mME$ID <- factor(new.mME$ID, levels = c('S1002_01_Non-responder', 'S1004_01_Non-responder', 'S1028_01_Non-responder', 'S1006_01_Responder',
                                         'S1013_01_Responder', 'S1014_01_Responder', 'S1020_01_Responder', 'S1037_01_Responder'))
## Relate modules to life stage  

#Prepare trait data. Data has to be numeric, so I will substitute time points/developmental stages for numeric values. 
#The "trait" we are considering here is lifestage. Make a dataframe that has a column for each lifestage name and a row for samples. 
#Populate a 1 for samples that match each lifestage and a 0 for samples not matching respective lifestages. 
#This process changes lifestages from a categorical variable into a binary variable. 
#This will allow for correlations between mean eigengenes and lifestage. 
metadata$num <- c("1")
MDP.meta = metadata %>%
  filter(Group == 'MDP')
allTraits <- as.data.frame(pivot_wider(MDP.meta, names_from = Responder, values_from = num, id_cols = Sample))
allTraits[is.na(allTraits)] <- c("0")
rownames(allTraits) <- allTraits$ID
datTraits <- allTraits[,c(-1)]
head(datTraits)
#define numbers of metabolites and samples and view 
nMetabolites = ncol(data_filt)
nSamples = nrow(data_filt)
nMetabolites
nSamples
# Correlations of traits with eigengenes
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
Colors=sub("ME","", names(MEs))
moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average")
pdf(file="ModuleTraitClusterTree.pdf", height=8, width=22)
plot(moduleTraitTree)
dev.off()

# Correlations of metabolites with eigengenes. Calculate correlations between ME's and groups 
moduleGeneCor=cor(MEs,data_filt)
data_filt[1:5, 1:5]
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);
head(moduleGenePvalue)
## Plot module-trait associations

#Represent module trait correlations as a heatmap.  

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
textMatrix = paste(signif(moduleTraitPvalue,1))
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
order <- c('Non-responder', 'Responder')
moduleTraitCor <- 
  moduleTraitCor[, c('Non-responder', 'Responder')]
datTraits <- 
  datTraits[, c('Non-responder', 'Responder')]

pdf(file="Module-trait-relationships_nolegend.pdf", height=10, width=5)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), 
               cex.lab.y= 1, cex.lab.x= 1, colors = blueWhiteRed(50),  setStdMargins = TRUE, 
               cex.text = 1, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()
ggsave('sample.png')

# 3. Plot module trait associations with a complex heatmap   
library(ComplexHeatmap)
library(dendsort)
#bold sig p-values
#dendrogram with WGCNA MEtree cut-off
#colored y-axis
#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)
#Make list of heatmap row colors
htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

#row_ha = rowAnnotation(ModuleSize = anno_text("11", "127", "8", "64", "85"), just = "left", 
#        location = unit(0.5, "npc"), show_name = TRUE)
# brown (11), grey (127), yellow (8), blue (64), turqoise (85) #figure out how to do row annotations to add sample sizes

time_order<-c("Non-responder", "Responder")

pdf(file = "Module-trait-relationship-heatmapadvanced.pdf", height = 11, width = 8)
ht=Heatmap(moduleTraitCor, name = "Eigengene", column_title = "Module-Group Eigengene Correlation", 
           col = blueWhiteRed(50), 
           row_names_side = "left", 
           row_dend_side = "left",
           width = unit(1.5, "in"), 
           height = unit(4.5, "in"), 
           #column_dend_reorder = TRUE, 
           #cluster_columns = col_dend,
           row_dend_reorder = FALSE,
           #column_split = 4, 
           row_split = 3, 
           #column_dend_height = unit(.5, "in"),
           cluster_rows = row_dend, 
           column_order = time_order, 
           row_gap = unit(1, "mm"), 
           border = TRUE,
          # cell_fun = function(j, i, x, y, w, h, col) {
          #   if(heatmappval[i, j] < 0.05) {
          #     grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
          #   }
          #   else {
          #     grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
          #   }},
           column_names_rot = 45,
           column_names_gp =  gpar(fontsize = 14, color='black', face='bold', border=FALSE),
           row_names_gp = gpar(fontsize = 14, color='black', face= 'bold', border = FALSE))
draw(ht)
dev.off()
# 4. Plot mean eigengene values for each module

#Load metadata and plot expression plots.  
# Number of each metabolites in each module 
row_mod_num <- data.frame(table(module_df$colors)) 

# Bring in additional metadata

MDP.meta$ID <- as.factor(MDP.meta$ID)
MDP.meta$Responder <- as.factor(MDP.meta$Responder)
MDP.meta$ID = paste(MDP.meta$Sample, MDP.meta$Responder, sep = '_')
mME_meta <- merge(mME, MDP.meta, by = "ID") %>%
  rename(Module = name)
mME_meta<-mME_meta%>%
  mutate(Cluster=case_when(Module %in% c('14', '1', '24', '18', '5', '22', '2', '4', '11', '12', '19', '3', '26') ~ 'Cluster1',
                           Module %in% c('10', '7', '9', '13', '17', '28', '23', '6') ~ 'Cluster2',
                           TRUE ~'Cluster3'))
mME_meta$Cluster<-as.factor(mME_meta$Cluster)

#Plot by time point with individual modules.  

mME_meta$Responder <- factor(mME_meta$Responder, levels = c("Non-responder", "Responder"))
mME_meta$Module <- factor(mME_meta$Module, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8",
                                                      '9',"10", "11", "12", "13", "14", "15", "16", "17", "18",
                                                      '19', "20", "21", "22", "23", "24", '25', '26', '27', '28'))
expression_plots<-mME_meta%>%
#  group_by(Module) %>%
  ggplot(aes(x=Responder, y=value, colour=Responder)) +
  facet_wrap(~ Module)+
  ylab("Module Expression") +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5), alpha = 0.7) +
  stat_summary(fun=mean, geom="line", aes(group=Responder, color = Responder), position = position_dodge(width = 0.5))  + 
  geom_point(pch = 21, position = position_dodge(width = 0.5)) +
  scale_fill_manual(name="Responder", values=c("#8C510A", "#DFC27D")) +
  scale_colour_manual(name="Responder", values=c("#8C510A", "#DFC27D")) + 
  xlab("") + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 11 , color = "black"),
        axis.title = element_text(size = 16, color = "black"), 
        axis.text.x = element_text(size=11, color="black", angle=45, hjust=1), 
        legend.title=element_blank(), 
        legend.text=element_text(color="black", size=12)); expression_plots

ggsave(filename="expression_eigengene.jpeg", plot=expression_plots, dpi=300, width=12, height=12, units="in")
head(mME_meta)
write.csv(mME_meta, "module_expression.csv")

#Add pattern names to each cluster. Cluster 1 = peaking in mid development; Cluster 2 = increasing; Cluster 3 = decreasing.  

mME_meta<-mME_meta%>%
  mutate(pattern=if_else(Cluster=="Cluster1", "Increasing", 
                         if_else(Cluster=="Cluster2", "Decreasing", 
                                 if_else(Cluster=="Cluster3", "NO_CHANGE", "NA"))))
mME_meta$pattern<-factor(mME_meta$pattern, levels=c("Decreasing", "Increasing", 'NO_CHANGE'))
#Display the number of metabolites and modules in each cluster.  
mME_meta%>%
  group_by(pattern)%>%
  summarise(module=length(unique(Module)))
mME_meta$Responder <- factor(mME_meta$Responder, levels = c("Non-responder", "Responder"))
mME_meta$Cluster <- factor(mME_meta$Cluster, levels = c("Cluster1", "Cluster2", "Cluster3"))
#mME_meta$lifestage <- factor(mME_meta$lifestage, levels = c("Egg (1 hpf)", "Embryo (5 hpf)", "Embryo (38 hpf)", "Embryo (65 hpf)", "Larvae (93 hpf) ", "Larvae (163 hpf)", "Larvae (183 hpf)", "Larvae (231 hpf) ", "Metamorphosed Polyp (183 hpf)", "Metamorphosed Polyp (231 hpf)", "Attached Recruit (183 hpf) ", "Attached Recruit (231 hpf)", "Attached Recruit (255 hpf)"))

expression_plots_cluster<-mME_meta%>%
  #  group_by(Module) %>%
  ggplot(aes(x=Responder, y=value)) +
  facet_grid(~pattern)+
  ylab("Eigengene Expression") +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_point(aes(fill=Responder, color=Responder), size=3, pch = 21, position = position_dodge(width = 0.5)) +
  geom_smooth(aes(group=1), se=TRUE, show.legend = NA, color="black")+
  scale_fill_manual(name="Responder", values=c("#8C510A", "#DFC27D")) +
  scale_color_manual(name="Responder", values=c("#8C510A", "#DFC27D")) +
  xlab("") + #Axis titles
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 11 , color = "black"),
        axis.title = element_text(size = 16, color = "black", face="bold"), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position="bottom", 
        legend.title=element_text(size=12, face="bold"), 
        legend.text=element_text(color="black", size=12), 
        strip.text.x = element_text(size = 14, color = "black", face="bold"), 
        strip.background = element_rect(color="white", fill="white")); expression_plots_cluster

ggsave(filename="expression_eigengene_cluster.jpeg", plot=expression_plots_cluster, dpi=300, width=10, height=5, units="in")

#############PATHWAY ANALYSIS
module_df <- module_df %>% rename(Module = colors)

# Select metabolites for each module 
cluster1<-c('14', '1', '24', '18', '5', '22', '2', '4', '11', '12', '19', '3', '26')
cluster2<-c('10', '7', '9', '13', '17', '28', '23', '6')


module_cluster1 <- module_df %>% filter(Module %in% cluster1)  %>% dplyr::select(Metabolite)
module_cluster2 <- module_df %>% filter(Module %in% cluster2)  %>% dplyr::select(Metabolite)
#module_cluster3 <- module_df %>% filter(Module %in% cluster3)  %>% dplyr::select(Metabolite)

#output to files for use
module_cluster1 <- module_cluster1%>%
  mutate(Metabolite = sub('^C_', '', Metabolite))
module_cluster2 <- module_cluster2%>%
  mutate(Metabolite = sub('^C_', '', Metabolite))

write.csv(module_cluster1, file="cluster1_metabolite_list.csv")
write.csv(module_cluster2, file="cluster2_metabolite_list.csv")


dim(module_cluster1)
dim(module_cluster2)


cluster1_kegg<-read.csv("cluster1_kegg_enrichment.csv")
cluster1_kegg$Cluster<-"Cluster1"

cluster2_kegg<-read.csv("cluster2_kegg_enrichment.csv")
cluster2_kegg$Cluster<-"Cluster2"

kegg_pathways<-rbind(cluster1_kegg, cluster2_kegg)
sig_pathways<-kegg_pathways%>%
  filter(FDR<0.06)%>%
  dplyr::select(X, total, hits, Raw.p, FDR, Cluster)

#plot significantly enriched terms by ontogenetic pattern 
Pathways.Dot.Plot <-  ggplot(sig_pathways, aes(x = Cluster, y = X)) + 
  geom_point(aes(colour=FDR, size=hits)) +
  #scale_color_continuous(name="FDR p-value", type = "gradient")+
  scale_color_distiller(type = "div",
                        direction = -1,
                        palette = "RdBu")+
  scale_size(range=c(4,11), limits=c(1,60))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
 # xlab("Cluster")+
 # ylab("KEGG Pathway")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 14, face="bold"),
                     strip.background = element_rect(fill="white"), 
                     axis.text.y = element_text(size = 14,  color = 'black'),
                     axis.text.x = element_text(size=14, face='bold', angle = 35, hjust=1, color = 'black'),
                    # axis.title.x = element_text(size = 14, face="bold"),
                    # axis.title.y = element_text(size = 14, face="bold"), 
                     legend.text = element_text(size=14), 
                     legend.title = element_text(size=14, face="bold"), 
                     legend.key.size = unit(0.2, 'cm'));Pathways.Dot.Plot

ggsave(filename="pathways_dot_plot.png", plot=Pathways.Dot.Plot, dpi=300, width=10, height=10, units="in")

#do same thing for compound class output
compound.class = read.csv('msea_pie_data._cluster1.csv')
head(compound.class)
c1.compound = compound.class%>%
  filter(Cluster =='Cluster1')
c2.compound = compound.class%>%
  filter(Cluster =='Cluster2')

ggplot(c2.compound, aes(x='', y= Hits, fill = Group))+
  geom_bar(width = 1, stat='identity', color = 'white')+
  coord_polar('y', start = 0)+
  #geom_text(aes(y= Hits, label = Hits), color = 'white')+
  scale_fill_manual(values = c2.compound$Color)+theme_void()
ggsave('look.png')

color = c('firebrick4', 'yellow', 'olivedrab', 'tan3', 'purple','thistle3', 'navyblue', 'gold', 'ivory4', 'black', 'plum2')
down$Group = factor(down$Group, levels = c("Lipids and lipid-like molecules", "Organic acids","Glycerolipids", 
                                       "Glycerophospholipids", "Fatty Acyls", "Organic oxygen compounds",
                                       "Polyketides" ,"Nucleic acids","Organoheterocyclic compounds", "Benzenoids", 
                                       "Sterol Lipids", 'Alkaloids', 'Carbohydrates', 'Organic nitrogen compounds'))

ggplot(up, aes(x='', y= Hits, fill = Group))+
  geom_bar(width = 1, stat='identity')+
  coord_polar('y', start = 0)+
  #geom_text(aes(y= Hits, label = Hits), color = 'white')+
  scale_fill_manual(values = color)+theme_void()
ggsave('look.png')

top = top[order(-top$Comp..1),]
top$X <- factor(top$X, levels = rev(top$X))
vip = ggplot(top, aes(x=X, y= Comp..1))+ 
  geom_point(size=5)+theme_minimal()+
  theme(panel.grid = element_blank(),  # Remove all grid lines
        panel.background = element_blank(),  # Set background to blank
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color='black', size=18),
        axis.text.y = element_text(color='black', size=18))+coord_flip()
ggsave('vip.png', dpi=300, width =15, height = 10)
