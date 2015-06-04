library(simpleaffy)
library(RColorBrewer)
library(affyPLM)
# Nuskaitomi .CEL failai is "celfiles" aplanko, pagal annotation.txt nurodyta tvarka.
celfiles <- read.affy(covdesc="annotation.txt", path="celfiles")
# Normalizuojami .CEL failai.
celfiles.gcrma <- gcrma(celfiles)
# Parenkama spalvu palete.
cols <- brewer.pal(8, "Set1")
# Nenormalizuotu intensyvumo reiksmiu grafikas.
png("All/bpunnorml.png", units="px", width=1600, height=1600, res=300)
par(mar = c(8, 2, 0.1, 0.001))
boxplot(celfiles, col=cols, las=2)
dev.off()
# Normalizuotu intensyvumo reiksmiu grafikas.
png("All/bpnorml.png", units="px", width=1600, height=1600, res=300)
par(mar = c(8, 2, 0.1, 0.001))
boxplot(celfiles.gcrma, col=cols, las=2)
dev.off()
# Tankio ir logaritminio intensyvumo histograma nenormalizuotiems duomenims.
png("All/histunorml.png", units="px", width=1600, height=1600, res=300)
hist(celfiles, col=cols)
dev.off()
# Tankio ir logaritminio intensyvumo histograma normalizuotiems duomenims.
png("All/histnorml.png", units="px", width=1600, height=1600, res=300)
hist(celfiles.gcrma, col=cols)
dev.off()
# probe-level metric skaiciavimai CEL failams.
celfiles.qc <- fitPLM(celfiles)
# affyPLM pateikia informatyvesnius grafikus
# RLE (Relative Log Expression) reiksmes turetu buti artimos nuliui.
png("All/RLE.png", units="px", width=1600, height=1600, res=300)
par(mar = c(8, 2.5, 3, 0.5))
RLE(celfiles.qc, las=2, main="RLE")
dev.off()
# NUSE (Normalised Unscaled Standard Errors).
# Medianos standartine paklaida turetu buti 1 daugumai genu.
png("All/NUSE.png", units="px", width=1600, height=1600, res=300)
par(mar = c(8, 2.5, 3, 0.5))
NUSE(celfiles.qc, las=2, main="NUSE")
dev.off()

celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=TRUE, remove.dupEntrez=TRUE)
# Kas isfiltruojama, jeigu nera Entrez gene ID
celfiles.filtered$filter.log

samples <- celfiles.gcrma$Factors
samples <- as.factor(samples)
design <- model.matrix(~0 + samples)
colnames(design) <- c("ACs1", "ACs2", "SCCs1", "SCCs2")

library(limma)
# Pritaikomas linijinis metodas genu raiskos rinkiniu filtravimui.
fit <- lmFit(exprs(celfiles.filtered$eset), design)
# Kontrasto matrica kiekvienos vezio grupes (ACs1, ACs2, SCCs1, SCCs2) palyginimui vienai su kita.
contrast.matrix <- makeContrasts(ACs1-ACs2, ACs1-SCCs1, ACs1-SCCs2, ACs2-SCCs1, ACs2-SCCs2, SCCs1-SCCs2, levels=design)
# Kontrasto matrica sujungiama su linijiniais metodais filtruotomis genu raiskomis.
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

library(hgu133plus2.db)
library(annotate)
require(ggplot2)
# coef=1 yra ACs1-ACs2, coef=2 yra ACs1-SCCs1, coef=3 yra ACs1-SCCs2
# coef=4 yra ACs2-SCCCs1, coef=5 yra ACs2-SCCs2, coef=6 yra SCCs1-SCCs2

for (j in 1:6){
probeset.list <- topTable(ebFit, coef=j, number=10000, lfc=2)
gene.symbols <- getSYMBOL(row.names(probeset.list), "hgu133plus2")
results <- cbind(probeset.list, gene.symbols)
gene_list <- topTable(ebFit, coef=j, number=1000000, sort.by="logFC")
gene.symbols <- getSYMBOL(row.names(gene_list), "hgu133plus2")
gene_list <- cbind(gene_list, gene.symbols)
par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
no_of_genes = dim(probeset.list)[1]
gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)

if (j == 1){write.csv(results, "All/ACs1vsACs2.csv", row.names=TRUE)
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + theme(legend.position = "none") + xlim(c(-5, 5)) + ylim(c(0, 7.5)) + xlab("log2 fold change") + ylab("-log10 p-value")
png("All/ACs1vsACs2.png", units="px", width=1600, height=1600, res=300)
plot(g)
dev.off()}

if (j == 2){write.csv(results, "All/ACs1vsSCCs1.csv", row.names=TRUE)
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + theme(legend.position = "none") + xlim(c(-5, 5)) + ylim(c(0, 7.5)) + xlab("log2 fold change") + ylab("-log10 p-value")
png("All/ACs1vsSCCs1.png", units="px", width=1600, height=1600, res=300)
plot(g)
dev.off()}

if (j == 3){write.csv(results, "All/ACs1vsSCCs2.csv", row.names=TRUE)
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + theme(legend.position = "none") + xlim(c(-5, 5)) + ylim(c(0, 7.5)) + xlab("log2 fold change") + ylab("-log10 p-value")
png("All/ACs1vsSCCs2.png", units="px", width=1600, height=1600, res=300)
plot(g)
dev.off()}

if (j == 4){write.csv(results, "All/ACs2vsSCCs1.csv", row.names=TRUE)
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + theme(legend.position = "none") + xlim(c(-5, 5)) + ylim(c(0, 7.5)) + xlab("log2 fold change") + ylab("-log10 p-value")
png("All/ACs2vsSCCs1.png", units="px", width=1600, height=1600, res=300)
plot(g)
dev.off()}

if (j == 5){write.csv(results, "All/ACs2vsSCCs2.csv", row.names=TRUE)
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + theme(legend.position = "none") + xlim(c(-5, 5)) + ylim(c(0, 7.5)) + xlab("log2 fold change") + ylab("-log10 p-value")
png("All/ACs2vsSCCs2.png", units="px", width=1600, height=1600, res=300)
plot(g)
dev.off()}

if (j == 6){write.csv(results, "All/SCCs1vsSCCs2.csv", row.names=TRUE)
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + theme(legend.position = "none") + xlim(c(-5, 5)) + ylim(c(0, 7.5)) + xlab("log2 fold change") + ylab("-log10 p-value")
png("All/SCCs1vsSCCs2.png", units="px", width=1600, height=1600, res=300)
plot(g)
dev.off()}
}