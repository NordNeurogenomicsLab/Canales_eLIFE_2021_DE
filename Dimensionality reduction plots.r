min.cpm.criteria <- 0.1 # really relaxed
test.data <- exp.data
#colnames(test.data) <- paste(sample.info$ExperimentalDesign, 1:102, sep = "_")
test.samples <- 1:nrow(sample.info)
min.cpm <- min.cpm.criteria
y <- DGEList(counts=test.data, group=ifelse(sample.info[,"Condition"]=="Saline",1,2))
keep <- rowSums(cpm(y)>min.cpm) >=2 #keeps only genes expressed in above min.cpm in at least 2 libraries in each group 
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y) # Normalizes for RNA composition
y <- estimateCommonDisp(y) # Estimates common dispersions. Calculates pseudo-counts, a type of normalized counts. Don't misteke them with 0+x type counts from other packages. Also "users are advised not to interpret the psuedo-counts as general-purpose normalized counts".
y <- estimateTagwiseDisp(y) #Estimates dispersions. Applicable only to experiments with single factor.
#Alternatively the two commands from above can be replaced with y <- estimateDisp(y)


#PCA plot

global.pseudo <- y$pseudo.counts
rownames(global.pseudo) <- rownames(y$counts) # this doesn't feel necessary because all(rownames(global.pseudo) == rownames(y$counts))
pca.results <- prcomp(scale(log(global.pseudo+1), center=T, scale=T)) # It is a good idea to scale your variables. Otherwise the magnitude to certain variables dominates the associations between the variables in the sample.

b <- data.frame(Samples=rownames(pca.results$rotation), 
                DPC=as.character(sample.info$DPC), 
                ExperimentalDesign=sample.info$ExperimentalDesign, 
                Condition=sample.info$Condition,
                sex.by.rna=sample.info$sex.by.rna,
                rRNA=sample.info$rRNA,
                lane=as.factor(sample.info$Lane),
                pca.results$rotation)
rownames(b) <- NULL

b$DPC <- ifelse(b$DPC == 19.5, "P0", b$DPC)

#PCA plot
PCA_plot <- 
  ggplot(b, aes(x=PC1, y=PC2, fill = DPC, colour=DPC, shape=Condition))+
  geom_point(size=3, alpha=0.6, aes(shape=Condition))+
  theme_bw()+
  labs(title = "PCA plot")+
  theme(plot.title = element_text(size = rel(1.5), hjust=0.5))

PCA_plot


#PCA plot lacking lane 12 and the inhibitor samples - slightly different version
PCA_plot_sex <- ggplot(b, aes(x=PC1, y=PC2, colour=sex.by.rna, shape=Condition))+
  geom_point(size=3, alpha=0.6)+
  theme_bw()+
  labs(title = "PCA plot")+
  theme(plot.title = element_text(size = rel(1.5), hjust=0.5))
  

PCA_plot_sex
