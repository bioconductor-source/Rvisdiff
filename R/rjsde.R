DEreport <- function(DE, counts = NULL, groups = NULL,
    cutoff = 0.05, normalized = NULL, directory = "DEreport"){

    if(!is.null(normalized)){
        counts <- normalized
        normalized <- TRUE
    }else{
        normalized <- FALSE
    }

    DE <- handleDEbyClass(DE, counts, cutoff, normalized, directory)

    if(identical(DE,FALSE)){
        return(invisible(NULL))
    }else{
        createReport(DE, counts, groups, cutoff, normalized, directory)
    }
}

createReport <- function(DE, counts, groups, cutoff, normalized, directory){
    create_directory(directory)

    groups <- checkGroups(DE,groups)

    if(!inherits(DE,'data.frame')){
        DE <- as.data.frame(DE)
    }

    for(p in c("pvalue","P.Value","PValue","padj","FDR","adj.P.Val")){
        if(p %in% colnames(DE)){
            DE[[p]][is.na(DE[[p]])] <- 1
        }
    }

    genes <- rownames(DE)
    if(sum(duplicated(genes))){
        warning("DE: some gene names are duplicated")
    }
    if(!"genes" %in% colnames(DE)){
        DE <- cbind(genes=genes,DE)
    }

    CPM <- NULL
    if(!is.null(counts)){
        if(!is.null(groups)){
            counts <- counts[,order(groups)]
            groups <- groups[order(groups)]
        }
        if(normalized){
            CPM <- counts[genes,]
        }else{
            CPM <- cpm(counts[genes,])
        }
    }

    createHTML(DE, CPM, cutoff, groups, normalized, directory)
}

checkGroups <- function(DE,groups){
    if(is.null(groups)){
        if(inherits(DE,"DESeqDataSet")){
            groups <- colData(DE)$conditions
        }else if(inherits(DE,"DGEList")){
            groups <- DE$samples$group
        }else if(inherits(DE,"DGEExact") || inherits(DE,"DGELRT")){
            if(!is.null(groups) &&
                !identical(setdiff(groups, DE$comparison), character(0))){
                groups <- NULL
                warning("groups: there are unmatching groups")
            }
        }else if(inherits(DE,"MArrayLM")){
            # TODO
        }
    }
    return(groups)
}

getColnamesJSON <- function(colnamesDE){
    pvalue <- "pvalue"
    if("P.Value" %in% colnamesDE)
        pvalue <- "P.Value"
    if("PValue" %in% colnamesDE)
        pvalue <- "PValue"

    padj <- "padj"
    if("FDR" %in% colnamesDE)
        padj <- "FDR"
    if("adj.P.Val" %in% colnamesDE)
        padj <- "adj.P.Val"

    log2FC <- "log2FoldChange"
    if("logFC" %in% colnamesDE)
        log2FC <- "logFC"

    expMean <- "baseMean"
    if("AveExpr" %in% colnamesDE)
        expMean <- "AveExpr"
    if("logCPM" %in% colnamesDE)
        expMean <- "logCPM"

    stat <- "stat"
    if("t" %in% colnamesDE)
        stat <- "t"
    if("LR" %in% colnamesDE)
        stat <- "LR"

    return(paste0('["', expMean, '","', log2FC, '","',
        pvalue, '","', padj, '","', stat, '"]'))
}

createHTML <- function(DE, CPM, cutoff, groups, normalized, directory){
    if(is.null(groups)){
        groups <- 'false'
    }else{
        groups <- paste0('["',paste0(groups,collapse='","'),'"]')
    }
    if(normalized){
        normalized <- 1
    }else{
        normalized <- 0
    }
    json <- paste0('{"DE":', tableJSON(DE), ',"cpms":', tableJSON(CPM),
        ',"names":', getColnamesJSON(colnames(DE)), ',"cutoff":', cutoff,
        ',"groups":', groups, ',"normalized":', normalized, '}')
    www <- wwwDirectory()
    html <- scan(file = paste0(www, "/template.html"), what = character(0),
        sep = "\n", quiet = TRUE)
    html <- sub("<!--json-->",
        paste0('<script type="application/json" id="data">',
        json, '</script>'), html)

    dir.create(paste0(directory,"/js"),FALSE)
    dir.create(paste0(directory,"/css"),FALSE)
    dir.create(paste0(directory,"/images"),FALSE)
    for(i in seq_len(nrow(dependencies))){
        file.copy(paste0(www,"/",dependencies[i,1]),
            paste0(directory,"/",dependencies[i,2]))
    }
    write(html,paste0(directory,"/index.html"))
    msg <- paste0("The report has been generated in the \"",
        normalizePath(directory),"\" path.")
    message(msg)
}

handleDEbyClass <- function(DE, counts, cutoff, normalized, directory){
    if(inherits(DE,"DESeqDataSet")){
        DE <- handleDESeqDataSet(DE,counts,cutoff,normalized,directory)
    }else if(inherits(DE,"DGEList")){
        DE <- handleDGEList(DE,counts,cutoff,normalized,directory)
    }else if(inherits(DE,"DGEExact") || inherits(DE,"DGELRT")){
        DE <- topTags(DE, n=nrow(DE), adjust.method="BH", sort.by="none")
    }else if(inherits(DE,"MArrayLM")){
        DE <- handleMArrayLM(DE,counts,cutoff,normalized,directory)
    }
    return(DE)
}

handleDESeqDataSet <- function(DE,counts,cutoff,normalized,directory){
        coldata <- colData(DE)
        datanames <- head(colnames(coldata),-1)
        keepdata <- vapply(datanames,function(x,data){
            vec <- data[[x]]
            return(length(levels(vec))>1 && length(levels(vec))<length(vec))
        }, logical(1), data = coldata)
        datanames <- datanames[keepdata]
        if(length(datanames)>1 ||
            length(levels(coldata[[datanames[1]]]))>2){
            create_directory(directory)
            nav <- character(0)
            for(d in datanames){
                cmb <- combn(levels(coldata[[d]]),2)
                for(i in seq_len(ncol(cmb))){
                    contrast <- c(d,cmb[1,i],cmb[2,i])
                    subsamples <- rownames(coldata)[
                        coldata[[d]] %in% contrast[2:3]]
                    subgroups <- coldata[subsamples,d]
                    subsamples <- subsamples[order(subgroups)]
                    subgroups <- subgroups[order(subgroups)]
                    subcounts <- NULL
                    if(!is.null(counts))
                        subcounts <- counts[,subsamples]
                    createReport(results(DE, independentFiltering = FALSE,
                        contrast = contrast), subcounts, subgroups, cutoff,
                        normalized, paste0(directory,"/",
                        paste0(contrast,collapse="_")))
                    nav <- c(nav,paste0(contrast,collapse="_"))
                }
            }
            metaIndex(nav,directory)
            return(FALSE)
        }
        return(results(DE,independentFiltering=FALSE))
}

handleDGEList <- function(DE,counts,cutoff,normalized,directory){
    datanames <- levels(DE$samples$group)
    if(length(datanames)>2){
        create_directory(directory)
        cmb <- combn(datanames,2)
        nav <- character(0)
        for(i in seq_len(ncol(cmb))){
            cmbs <- c(cmb[1,i],cmb[2,i])
            subsamples <- rownames(DE$samples)[DE$samples$group %in% cmbs]
            subgroups <- DE$samples[subsamples,'group']
            subsamples <- subsamples[order(subgroups)]
            subgroups <- subgroups[order(subgroups)]
            subcounts <- NULL
            if(!is.null(counts))
                subcounts <- counts[,subsamples]
            edger <- exactTest(DE,pair=cmbs)
            edger <- topTags(edger, n=nrow(DE), 
                adjust.method="BH", sort.by="none")
            createReport(edger, subcounts, subgroups, cutoff, normalized,
                paste0(directory,"/",paste0(cmbs,collapse="_")))
            nav <- c(nav,paste0(cmbs,collapse="_"))
        }
        metaIndex(nav,directory)
        return(FALSE)
    }
    DE <- exactTest(DE,pair=c(1,2))
    DE <- topTags(DE, n=nrow(DE), adjust.method="BH", sort.by="none")
    return(DE)
}

handleMArrayLM <- function(DE,counts,cutoff,normalized,directory){
        datanames <- colnames(DE$design)
        if(length(datanames)>1){
            create_directory(directory)
            nav <- character(0)
            if(!is.null(DE$contrasts)){
                if(!is.null(counts))
                    samples <- colnames(counts)
                else
                    samples <- as.character(seq_len(nrow(DE$design)))
                groups <- character(length(samples))
                datanames <- rownames(DE$contrasts)[
                    as.logical(DE$contrasts[,1])]
                for(i in seq_along(samples)){
                    for(j in datanames)
                        if(DE$design[i,j]==1)
                            groups[i] <- j
                }
                names(groups) <- samples
                for(i in seq_len(ncol(DE$contrasts))){
                    colname <- colnames(DE$contrasts)[i]
                    cmb <- rownames(DE$contrasts)[DE$contrasts[,i]!=0]
                    subsamples <- samples[seq_len(nrow(DE$design))[
                        as.logical(DE$design[,cmb[1]]+DE$design[,cmb[2]])]]
                    subgroups <- groups[subsamples]
                    subsamples <- subsamples[order(subgroups)]
                    nav <- topTableReport(DE, nav, colname, i, counts,
                        subgroups, cutoff, normalized, directory)
                }
            }else{
                for(i in seq_along(datanames)){
                    colname <- datanames[i]
                    subgroups <- c("other",colname)
                    subgroups <- subgroups[(DE$design[,colname]+1)]
                    nav <- topTableReport(DE, nav, colname, i, counts,
                        subgroups, cutoff, normalized, directory)
                }
            }
            metaIndex(nav,directory)
            return(FALSE)
        }
        return(topTable(DE, coef = 2, number = nrow(DE),
            sort.by = "none", adjust.method = "BH"))
}

topTableReport <- function(DE, nav, colname, i, counts, subgroups, cutoff,
    normalized, directory){
    subcounts <- NULL
    if(!is.null(counts))
        subcounts <- counts[,order(subgroups)]
    subgroups <- subgroups[order(subgroups)]
    createReport(topTable(DE, coef = i, number = nrow(DE),
        sort.by = "none", adjust.method = "BH"), subcounts,
        subgroups, cutoff, normalized,
        paste0(directory,"/",colname))
    return(c(nav,colname))
}

metaIndex <- function(nav,directory){
    www <- wwwDirectory()
    html <- scan(file = paste0(www, "/meta.html"), what = character(0),
        sep = "\n", quiet = TRUE)

    nav <- paste0("<li><a href=\"", nav, "/index.html\">",
        nav, "</a></li>", collapse="")
    html <- sub("<!--nav-->",paste0("<ul>",nav,"</ul>"),html)

    dir.create(paste0(directory,"/css"),FALSE)
    for(i in c(5,6)){
        file.copy(paste0(www,"/",dependencies[i,1]),
            paste0(directory,"/",dependencies[i,2]))
    }
    write(html,paste0(directory,"/index.html"))
    msg <- paste0("The index has been generated in the \"",
        normalizePath(directory),"\" path.")
    message(msg)
}

dependencies <- data.frame(name=c("jquery-3.3.1.min.js",
    "d3.v3.min.js",
    "datatables.min.js",
    "graphs.js",
    "responsee.css",
    "graphs.css",
    "datatables.css",
    "select.datatables.css",
    "sorting.png",
    "sorting_asc.png",
    "sorting_asc_disabled.png",
    "sorting_desc.png",
    "sorting_desc_disabled.png"),
    type=c("js","js","js","js","css","css","css","css",
    "images","images","images","images","images"))


tableJSON <- function(x){
    if(is.null(x)){
        return('false')
    }
    for(col in colnames(x)){
        if(!inherits(x[,col],"numeric")){
            x[,col] <- paste0('"',x[,col],'"')
        }
    }
    colNames <- paste0('[',paste0('"',colnames(x),'"',collapse=','),']')
    if(inherits(x,'data.frame')){
        aux <- vapply(x,format,character(nrow(x)),trim=TRUE,justify="none")
    }else{
        aux <- format(x,trim=TRUE,justify="none",scientific=FALSE,digits=3)
    }
    aux[aux=="NA"] <- "null"
    aux <- apply(aux,1,function(x) paste0('[',paste0(x,collapse=','),']'))
    aux <- paste0(c(colNames,aux), collapse = ",")
    json <- paste0("[", aux, "]", collapse = "")
    return(json)
}

create_directory <- function(directory){
    if(file.exists(directory)){
        unlink(directory, recursive = TRUE)
    }
    dir.create(directory)
}

wwwDirectory <- function(){
    path <- system.file("www",package="Rvisdiff")
    return(path)
}

