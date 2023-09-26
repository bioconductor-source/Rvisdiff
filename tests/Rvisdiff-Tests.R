library(Rvisdiff)

### Prapare input data

library("airway")
data("airway")
se <- airway
se$dex <- relevel(se$dex, ref="untrt")
countdata <- assay(se)
coldata <- colData(se)

### Produce interactive DE reports from DESeq2 results

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- DESeq(dds)
DEreport(dds, countdata)

### Produce interactive DE reports from edgeR results

library("edgeR")
design <- model.matrix(~ cell + dex, data = coldata)
dl <- DGEList(counts = countdata, group = coldata$dex)
dl <- calcNormFactors(dl)
dl <- estimateDisp(dl, design=design)
DEreport(dl, countdata)
de <- exactTest(dl,pair=1:2)
DEreport(de, countdata, coldata$dex)

### Produce interactive DE reports from limma results

library("limma")
design <- model.matrix(~ 0 + dex + cell, data = coldata)
contr <- makeContrasts(dextrt - dexuntrt,levels=colnames(design))
limmaexprs <- voom(countdata, design)
fit <- lmFit(limmaexprs, design)
fit <- contrasts.fit(fit, contrasts=contr)
fit <- eBayes(fit)
DEreport(fit, countdata)

### Produce interactive DE reports from Differential test results

untrt <- countdata[,coldata$dex=="untrt"]
trt <- countdata[,coldata$dex=="trt"]
library(matrixTests)
wilcox <- col_wilcoxon_twosample(t(untrt), t(trt))
stat <- wilcox$statistic
p <- wilcox$pvalue
log2FoldChange <- log2(rowMeans(trt)+1) - log2(rowMeans(untrt)+1)
wilcox <- data.frame(genes = rownames(countdata), statistic = stat,
    pValue = round(p, 6), pAdjust = p.adjust(wilcox[,2], method = "BH"),
    expMean = rowMeans(countdata), log2FC = log2FoldChange)
normalized <- edgeR::cpm(countdata)
normalized <- as.data.frame(normalized)
normalized$genes <- rownames(normalized)
DEreport(wilcox, NULL, coldata$dex, normalized = normalized, genes="genes",
pvalue="pValue", padj = "pAdjust", stat = "statistic",
baseMean="expMean", log2FoldChange="log2FC")

### Missing columns tests
colnames(wilcox) <- c("genes","missing1","missing2","missing3","missing4","missing5")
DEreport(wilcox, NULL, coldata$dex, normalized = normalized, genes="genes")

