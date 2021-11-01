############################################################################################
######  STAMP Revisions Skin Control data analysis by Ximo Pechuan i Jorge 29/07/2021
############################################################################################

#################### Set-up
source("/Users/pechuanj/Documents/Workstation/stamp/STAMP_ngs_analyses/src/BulkRNAseq_AnalysisFunctions.R")
path.to.data = "/Users/pechuanj/Documents/Workstation/stamp/NGS3947/"
setwd(path.to.data)

################### Loading and prepare data
# Read sample annotations
pheno = read.csv(paste(path.to.data,"data/","NGS3947_pheno.csv",sep=""),row.names = 1)
# Separate sample label into its components
pheno = pheno %>%  separate(SAMPLE_LABEL,"_", into =c("ID1","Tissue","Day"))
pheno$Complete = paste(pheno$Tissue,pheno$Day,sep="_")
# Load eset 
es.all = readRDS("data/NGS3947_ESET.Rdata")
# DGEList 
y.all = DGEList(
  exprs(es.all),genes=fData(es.all))
y.all$samples = cbind(y.all$samples, pheno)
y.all = calcNormFactors(y.all)

############################ Filtere lowly expressed genes
keep.exprs = filterByExpr(y.all, group=y.all$samples$Complete)
yf =  y.all[keep.exprs,, keep.lib.sizes=FALSE]
yf = calcNormFactors(yf)

##########################################PCA plot
ExprMat = cpm(yf,log=T)
var_genes = apply(ExprMat, 1, var)
var_gens = var_genes[order(var_genes,decreasing = T)]
nvar = 500 # How many variable genes?
topgenes = var_genes[1:nvar]
topdat = ExprMat[names(topgenes),]

####################################### PCA
res.pca = PCA(t(topdat), scale.unit = TRUE, ncp = 10, graph = F)

pdf("PCA_NGS3947.pdf")
fviz_pca_ind(res.pca,
             axes = c(1,2),
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.character(pheno$Complete), # color by groups
             palette = "viridis", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black", repel = TRUE,
             legend.title = "",
             ellipse.level  = 0.90,
             pointshape = 21
) 
dev.off()

# Correlation
vint1 = p1$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vint2 = p2$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vint3 = p3$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vintmes = unique(c(as.character(vint1$name),as.character(vint2$name),as.character(vint3$name)))
correlacions = res.pca$var$cor[vintmes,1:3]

# Plot
pdf("Correlations.pdf")
ggcorrplot(t(correlacions), tl.cex  = 7)
dev.off()
################################################### PCA 3D #######
pca = prcomp(t(topdat), scale.=TRUE)
gr = factor(as.character(pheno$Complete))
library("car")
library("rgl")
scatter3d(pca$x[,1],pca$x[,2],pca$x[,3],groups = gr,
          surface=FALSE, grid = FALSE,
          ellipsoid = T,
          # surface.col = "viridis",
          xlab = "PC1", ylab = "PC2",
          zlab = "PC3",
          
          axis.col = c("black", "black", "black")
)
#rgl.postscript("Bulk2_3D.eps",fmt="eps")

############################ Progeny pathways
annotations = pheno[,c("Tissue","Day")]

# progeny
pathways = progeny(ExprMat,scale = T, organism = "Mouse", top = 100,
                   perm = 1, verose = FALSE)
myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
progenyHM = pheatmap::pheatmap(t(pathways),fontsize=14, show_rownames = T,
                               color=myColor, main = "PROGENy  Pathways", 
                               angle_col = 45, treeheight_col = 1,
                               show_colnames = F,
                               border_color = NA,
                               annotation_col = annotations
                               #annotation_colors = ann_colors
                               )

save_pheatmap_pdf(progenyHM, filename = paste("NGS3947","HeatmapProgeny.pdf",sep=""))

############################ Score/PCA by Inf signature both Alpha and Gamma
jeremie_inf = c("Rsad2",
                "Mx1",
                "Mx2",
                "Irf7",
                "Irf9",
                "Cxcl9",
                "Cxcl10",
                "Mb21d1",
                "Oasl2",
                "Oas1b",
                "Ifit1",
                "Ifit2",
                "Ifit3",
                "Ifitm1",
                "Ifitm2")
jeremie_dat = ExprMat[(rownames(ExprMat) %in% jeremie_inf),]
# PCA
res.jere = PCA(t(jeremie_dat), scale.unit = TRUE, ncp = 10, graph = F)
pdf("Jeremie_inf_genes.pdf")
fviz_pca_ind(res.jere,
             axes = c(1,2),
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.character(pheno$Complete), # color by groups
             palette = "viridis", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black", repel = TRUE,
             legend.title = "",
             ellipse.level  = 0.90,
             pointshape = 21
) 
dev.off()
# Violoin plot
#Score
pcSig = gsScore(jeremie_dat)
yf$samples$pcSig = pcSig
pdf("Jermie_Ifn_gamma_Score.pdf")
  ggviolin(yf$samples, "Complete", "pcSig", fill = "Complete")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# Pheatmap
jeremie_hm = pheatmap::pheatmap(jeremie_dat,fontsize=14, show_rownames = T,
                               color=myColor, main = "Jeremie INF Signature", 
                               angle_col = 45, treeheight_col = 1,
                               show_colnames = F,
                               border_color = NA,
                               annotation_col = annotations,
                               cluster_rows = T,
                               cluster_cols = F,
                               scale = "row"
                               #annotation_colors = ann_colors
)

save_pheatmap_pdf(jeremie_hm, filename = paste("NGS3947","INF_jeremie.pdf",sep=""))


###################################### DE and pathway analysis
# Tumor vs Pore | Day
# Skin vs Pore Day_n (wound healing over time)
# Pore(t) -> Also Through intersection method
# Tumor(t) -> Also Through intersection method
#################################################################

############################## Tumor vs Pore | Day
setwd(path.to.data)
Treatment = factor(pheno$Complete)
design = model.matrix(~0+Treatment)
colnames(design) = levels(Treatment)
# Run voom 
vmf=voomWithQualityWeights(yf, design, plot=T)
# Run Limma
fit=lmFit(vmf, design)
# Find differentially expressed genes
contrast.matrix = makeContrasts("Day1" = Tumor_Day1 - Pore_Day1,
                                "Day3" = Tumor_Day3 - Pore_Day3,
                                "Day8" = Tumor_Day8 - Pore_Day8,
                                 levels=design)

fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
dt = decideTests(fit2,lfc = 0,p.value = 0.05)
summary(decideTests(fit2))
# input ready
species = "Mus musculus"
logfold = 0
Contrastes = colnames(decideTests(fit2))[colSums(abs(decideTests(fit2))) > 0]
Prefix = "Each_day"
# all plots
BulkPlots(ExprMat,Prefix,Contrastes,annotations,Treatment,species,logfold)

########################### Skin vs Pore Day_n (wound healing over time)
setwd(path.to.data)
Prefix = "Skin_vs_pore"
# Find differentially expressed genes
contrast.matrix = makeContrasts("Day1" = Pore_Day1  - Skin_Day8,
                                "Day3" = Pore_Day3 - Skin_Day8,
                                "Day8" = Pore_Day8 - Skin_Day8,
                                levels=design)

fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
dt = decideTests(fit2,lfc = 0,p.value = 0.05)
summary(decideTests(fit2))

# input ready
species = "Mus musculus"
logfold = 0
Contrastes = colnames(decideTests(fit2))[colSums(abs(decideTests(fit2))) > 0]
Prefix = "Skin_Pores_days"
# all plots
BulkPlots(ExprMat,Prefix,Contrastes,annotations,Treatment,species,logfold)

# Venn diagram
pair = abs(dt)
interPair= rownames(pair[rowSums(pair) == 2,])
write.csv(interPair,"Genes_significant_in_both_contrasts.csv")
pdf(paste(Prefix,"SignatureVenn.pdf",sep="_"))
vennDiagram(dt, circle.col=ImmunePhenotypes,cex = 0.6)
dev.off()

# input ready
species = "Mus musculus"
logfold = 0
Contrastes = colnames(decideTests(fit2))[colSums(abs(decideTests(fit2))) > 0]

# all plots
BulkPlots(ExprMat,Prefix,Contrastes,annotations,Treatment,species,logfold)

#################################  Time Series Spline Model
# Remove skin as there is no timeseries
pheno = pheno %>% filter(!(Tissue == "Skin"))
pheno$Day = fct_recode(pheno$Day,"1"="Day1","3"="Day3","8"="Day8")
# remove pheno
es.filt = es.all[,rownames(pheno)]
# DGEList 
y.all = DGEList(
  exprs(es.filt),genes=fData(es.filt))
y.all$samples = cbind(y.all$samples, pheno)
y.all = calcNormFactors(y.all)

############################ Filtere lowly expressed genes
lowExprs = rowSums(edgeR::cpm(es.filt) >= 0.2) < 5
yf =  y.all[!(lowExprs),, keep.lib.sizes=FALSE]
yf = calcNormFactors(yf)
# Factors
Tissue = factor(pheno$Tissue)
Time = factor(pheno$Day)

######################### Time Spline
# Fit a five degree polynomial spline
X = ns(as.numeric(as.character(Time)), df=2)
# Design matrix
desing_mat = model.matrix(~Tissue*X)
# Run voom 
vmf=voomWithQualityWeights(yf, desing_mat, plot=T)
# Run Limma
fit=lmFit(vmf, desing_mat)
fit2 = eBayes(fit)
plotSA(fit2, main="Final model: Mean-variance trend")
# Interaction genes
genes = topTable(fit2, coef=5:6,n=Inf,adjust.method="fdr",sort.by = "F")
genes$gene = rownames(genes)
genes = genes %>% filter(adj.P.Val<0.05)
dir.create("Tumor_Pore_Interaction_SplineMethod")
setwd("Tumor_Pore_Interaction_SplineMethod")
write.csv(genes, 'Tumor_Pore_interaction_genes.csv')
ExprMat = cpm(yf,log=T)
# Plot interacting genes
for (i in 1:length(genes$gene)) {
  gensmbl= genes$gene[i]
  filename = paste(gensmbl,"_","Interaction.pdf",sep="")
  pdf(filename)
  print(
    interaction.plot(x.factor = Time,    # variable to plot on x-axis
                     trace.factor = Tissue, # variable to specify "traces"; here, lines
                     response = ExprMat[gensmbl,],    # variable to plot on y-axis
                     fun = mean,  # summary statistic to be plotted for response variable
                     type = "l",     # type of plot, here "l" for lines
                     ylab = "Log(cpm)",
                     xlab = "Time (Hours)",
                     col = c("orange4", "red4","blue4"),
                     lty = 1,  # line type
                     lwd = 2,  # line widthline_reg_1953_ko.txt
                     trace.label = "Tissue",  # label for legend
                     xpd = T,
                     main=gensmbl) #,  # 'clip' legend at border
  )
  dev.off()
  
}

