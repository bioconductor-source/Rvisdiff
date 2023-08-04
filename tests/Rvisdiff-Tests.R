library(Rvisdiff)

### Produce interactive DE reports

library("airway")
data("airway")
se <- airway
se$dex <- relevel(se$dex, ref="untrt")
countdata <- assay(se)

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- DESeq(dds)
dr <- results(dds, independentFiltering = FALSE)

DEreport(dr, countdata, se$dex)

