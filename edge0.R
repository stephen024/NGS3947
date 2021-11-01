#########################################################################
######  STAMP Day 13 bulk RNAseq data analysis by Ximo Pechuan i Jorge 
#########################################################################

# Set-up
source("iSTAMP_functions.R")
#setwd("/Users/pechuanj/Documents/Workstation/stamp/STAMP_ngs_analyses/data/")

# Read sample annotations
pheno = read.csv('data/NGS3947_pheno.csv',row.names = 1)
count_m=read.csv('data/NGS3947_counts.csv',row.names = 1)
es.all=readRDS('data/NGS3947_ESET.Rdata')
label_in=pheno$SAMPLE_LABEL
sk_sam_index=grep("Skin",label_in)
po_sam_index=grep("Pore",label_in)

sk_po_pheno=pheno[c(sk_sam_index,po_sam_index),]
sk_po_count=count_m[,c(sk_sam_index,po_sam_index)]

f1=rep("Skin",times=length(sk_sam_index))
f2=rep("Pore",times=length(po_sam_index))

group=factor(c(f1,f2),levels = c("Skin","Pore"))
deg_obj=DGEList(counts = sk_po_count,group = group,genes=fData(es.all))
###. create edgeR object ##

lcpm=cpm(deg_obj,log = TRUE)

keep.exprs=filterByExpr(deg_obj,group=group)
deg_obj <- deg_obj[keep.exprs,, keep.lib.sizes=FALSE]
dim(deg_obj)


deg_obj <- calcNormFactors(deg_obj, method = "TMM")
plotMDS(lcpm,labels = group)

design<-model.matrix(~0+group)
colnames(design)=levels(group)
deg_obj= estimateDisp(deg_obj,design,robust = TRUE)
plotBCV(deg_obj)

con=makeContrasts(Skin-Pore,levels = design)
fit=glmQLFit(deg_obj,design,robust = TRUE)
head(fit$coefficients)
qlf=glmQLFTest( fit,contrast=con)

summary(decideTests(qlf))
plotMD(qlf)


es.all=readRDS('data/NGS3947_ESET.Rdata')
#a=fData(es.all)
#b=exprs(es.all)


y.all<-DGEList(exprs(es.all),genes=fData(es.all))

# y.all$samples
y.all=calcNormFactors(y.all)


# Filter bad quality samples
remove = c("SAM24376843","SAM24376842","SAM24376828",
           "SAM24376836","SAM24376835","SAM24376852")
pheno = pheno[!(rownames(pheno) %in% remove),]

pheno = pheno %>% filter(Immune %in% c("Excluded","Desert","Inflamed"))
# Keep only the immune phenotypes
pheno$tdTomato = log(pheno$tdTomato+1)
pheno$GFP = log(pheno$GFP+1)

# Load data
load('data/NGS3947_ESET.Rdata')
keep = rowSums(edgeR::cpm(es.all) >= 0.2) > 6
es.all = es.all[,rownames(pheno)]
es.filt = es.all[keep,]

# Change to symbol
annotations = es.filt@featureData@data
rownames(es.filt) = uniquifyFeatureNames(ID= rownames(es.filt),names=annotations[rownames(es.filt), "symbol"])
data = exprs(es.filt)

# Save the csv
write.csv(data,"bulk3_OnlyImmuneFilt_data.csv")
write.csv(pheno,"pheno_OnlyImmuneFilt_bulk3.csv")

#Prepare for PCA plot
gr <- factor(pheno$Immune)
data = exprs(es.filt)
logdat = log(data+1)
var_genes = apply(logdat, 1, var)
var_genes = var_genes[order(var_genes,decreasing = T)]
plot((var_genes[1:5000])) # By projection score
topgenes = var_genes[1:5000]
topdat = logdat[names(topgenes),]

# PCAw
res.pca = PCA(t(topdat), scale.unit = TRUE, ncp = 10, graph = F)

scatter3d(pca$x[,1],pca$x[,2],pca$x[,3],groups = gr,
          surface=FALSE, grid = FALSE, ellipsoid = F,
          surface.col = ImmunePhenotypes,
          xlab = "PC1", ylab = "PC2",
          zlab = "PC3",
          axis.col = c("black", "black", "black")
)
#rgl.postscript("Bulk2_3D.eps",fmt="eps")
rgl.postscript("plot.pdf",fmt="pdf")
rgl.snapshot(filename = "D13_only_Immune.png")

pdf("PC1_PC2_immune.pdf")
fviz_pca_ind(res.pca,
             axes = c(1,2),
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pheno$Immune, # color by groups
             palette = ImmunePhenotypes, 
             addEllipses = T, label = "var",
             col.var = "black", repel = TRUE,
             legend.title = "Immune Phenotype",
             title = "Day 13", ellipse.level  = 0.8
) 
dev.off()


# Contribution plots
p1 = fviz_contrib(res.pca, choice = "var", axes = 1, top = 20) + theme(axis.text=element_text(size=6))
p2 = fviz_contrib(res.pca, choice = "var", axes = 2, top = 20) + theme(axis.text=element_text(size=6))
p3 = fviz_contrib(res.pca, choice = "var", axes = 3, top = 20) + theme(axis.text=element_text(size=6))
p4 = fviz_contrib(res.pca, choice = "var", axes = 4, top = 20) + theme(axis.text=element_text(size=6))
pdf("Contributions_plot.pdf")
plot_grid(p1, p2, p3,p4,labels = c('', '',"",""))
dev.off()

# Correlation
vint1 = p1$data %>% arrange(-contrib) %>% top_n(contrib,n=20)
vint2 = p2$data %>% arrange(-contrib) %>% top_n(contrib,n=20)
vint3 = p3$data %>% arrange(-contrib) %>% top_n(contrib,n=20)
vint4 = p4$data %>% arrange(-contrib) %>% top_n(contrib,n=20)
vintmes = unique(c(as.character(vint1$name),as.character(vint2$name),
                   as.character(vint3$name),as.character(vint4$name)))
correlacions = res.pca$var$cor[vintmes,1:4]

# Plot
pdf("Correlations.pdf")
ggcorrplot(t(correlacions), tl.cex  = 10)
dev.off()

# Differential Expression Contrasts
Treatment = factor(pheno$Immune)
GFP = pheno$GFP
design = model.matrix(~0+Treatment)
colnames(design) <- c("Desert","Excluded","Inflamed")

vmf=voomWithQualityWeights(es.filt, design, plot=TRUE)


## Run Limma
fit=lmFit(vmf, design)

contrast.matrix = makeContrasts("ExcludedSignature" = Excluded -(Desert+Inflamed)/2,
                                 "InflamedSignature" = Inflamed -(Excluded+Desert)/2,
                                 "DesertesertSignature" = Desert -(Excluded+Inflamed)/2,
                                 "ExcludedvsInflamed" = Excluded - Inflamed,
                                 "DesertvsInflamed" = Desert - Inflamed,
                                "DesertvsExcluded" = Desert - Excluded,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
dt = decideTests(fit2,lfc = 0,p.value = 0.05)

# input ready
ExprMat = vmf$E
pheno$Immune = as.character(pheno$Immune)
annotations = pheno[,c("Immune","tdTomato","GFP")]
Prefix = "Bulk3_D13_onlyImmune"
species = "Mus musculus"
logfold = 0


summary(decideTests(fit2))
Contrastes = colnames(decideTests(fit2))

### Venn
pdf(paste(Prefix,"SignatureVenn.pdf",sep="_"))
vennDiagram(dt[,1:3], circle.col=ImmunePhenotypes,cex = 0.6)
dev.off()
pdf(paste(Prefix,"PairwiseVenn.pdf",sep="_"))
vennDiagram(dt[,c(5,6)], circle.col=ImmunePhenotypes,cex = 0.6)
dev.off()

## 
pair = abs(dt[,5:6])
interPair= rownames(pair[rowSums(pair) == 2,])
#
onlyimmune = annotations %>% filter(Immune %in% c("Desert","Excluded","Inflamed"))
onlyimmune$Immune = factor(onlyimmune$Immune)
signatura = interPair

my_heatmap = pheatmap::pheatmap(ExprMat[signatura,rownames(onlyimmune)], 
                                border_color = T,scale="row",fontsize_row=3 ,
                                show_colnames=F,
                                annotation_col = onlyimmune,annotation_colors = ann_colors)

save_pheatmap_pdf(my_heatmap, filename = paste(Prefix,"HeatmapPair.pdf",sep=""))
write.csv(signatura,"Pair_Contrasts_Discriminating_signature.csv")

## 
pair = abs(dt[,1:3])
interPair= rownames(pair[rowSums(pair) == 3,])
#
onlyimmune = annotations %>% filter(Immune %in% c("Desert","Excluded","Inflamed"))
onlyimmune$Immune = factor(onlyimmune$Immune)
signatura = interPair

my_heatmap = pheatmap::pheatmap(ExprMat[signatura,rownames(onlyimmune)], 
                                border_color = T,scale="row",fontsize_row=3 ,
                                show_colnames=F,
                                annotation_col = onlyimmune,annotation_colors = ann_colors)

save_pheatmap_pdf(my_heatmap, filename = paste(Prefix,"HeatmapSig.pdf",sep=""))
write.csv(signatura,"Sig_Contrasts_Discriminating_signature.csv")

# progeny
pathways <- progeny(ExprMat,scale = T, organism = "Mouse", top = 100,
                    perm = 1, verbose = FALSE)
myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
progenyHM = pheatmap::pheatmap(t(pathways),fontsize=14, show_rownames = T,
         color=myColor, main = "PROGENy", angle_col = 45, treeheight_col = 1,
         show_colnames = F,
         border_color = NA,annotation_col = onlyimmune,annotation_colors = ann_colors)

save_pheatmap_pdf(progenyHM, filename = paste(Prefix,"HeatmapProgeny.pdf",sep=""))
# all plots

BulkPlots(ExprMat,Prefix,Contrastes,annotations,Treatment,species,logfold)




