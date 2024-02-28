#### functions organoids paper

library(clusterProfiler)
library(org.Hs.eg.db)
go_analysis <- function (query, uni, name, organism="mouse", path) {
  if (organism=="mouse") {
    db = org.Mm.eg.db
  }
  else if (organism=="human"){
    db = org.Hs.eg.db
  }
  uni <- bitr(uni, fromType = "SYMBOL", toType = c("ENTREZID", "GO"),OrgDb = db, drop=F)
  uni <- uni$ENTREZID
  uni <- uni[!is.na(uni)]
  uni <- unique(uni)
  cat(name, "\n")
  g <- query
  g <- bitr(g, fromType = "SYMBOL", toType = c("ENTREZID", "GO"),OrgDb = db, drop=F)
  g <- g$ENTREZID
  g <- g[!is.na(g)]
  g <- unique(g)
  for (l in c("BP", "CC", "MF")) {
    ego <- enrichGO(gene          = g,
                    universe      = uni,
                    OrgDb         = db,
                    ont           = l,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE,
                    minGSSize = 3)
    if (!is.null(ego)) {
      if (! nrow(ego) == 0) {
        pdf(file.path(path, paste0("GOenrich_nodes_", name, "_", l, ".pdf")), 12, 12)
        print(dotplot(ego))
        dev.off() }
      if (l == "BP"){
        return(ego)
      }
      else{ cat(paste("no go enrichment", l, "\n"))}
      return("NA")
    }
    
    else{
      cat(paste("no go enrichment", l, "\n"))
      return("NA")
    }
  }
}

get_sample_condition_filter <- function(sd_tab, cent=F, batch=NULL,drop=F) {
  if (cent) {
    sd_tab <- cbind(sd_tab, batch=batch)
    f <- ! batch %in% c("810", "373")
  }
  else{
    f <- ! sd_tab[,"batch"] %in% c("810", "373")}
  if(drop==T){
    sd_tab <- sd_tab[f,]}
  sample <- sd_tab[,"batch"]
  sample[which(sample %in% c("810", "716", "717"))] <- "condition1"
  sample[which(sample == "807")] <- "condition2"
  sample[which(sample == "857")] <- "condition3"
  sample[which(sample %in% c("808", "809"))] <- "condition4"
  sample[which(sample %in% c("855", "856"))] <- "condition5"
  sample[which(sample %in% c("858", "859"))] <- "condition6"
  sample[which(sample %in% c("373", "375"))] <- "condition7"
  sample[which(sample %in% "374")] <- "condition8"
  sample[which(sample %in% "861")] <- "condition9"
  sample[which(sample %in% "RPE")] <- "condition10"
  sd_tab[,"batch"] <- sample
  sd_tab
}

plotexpression <- function (x, y, g, n, col = NULL, name = NULL, cluster = FALSE,
                            alpha = 0.5, types = NULL, cex = 3, ylim = NULL, leg = T,
                            lwd = 2.5, legendpos="topleft", samp=F, samppart=NULL, sampcol=NULL) {
  library(RColorBrewer)
  
  if (length(g) <= 8) {
    set1 <- brewer.pal(length(g), "Set1")
    if (length(g) >= 6) {
      set2 <- brewer.pal(length(g), "Set2")
      set1[6] <- set2[6]}
  }
  else {
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set1 <- sample(col_vector, length(g))
  }
  
  for ( i in 1:length(g)){
    cl <- unique(y[n])
    set.seed(111111)
    if (is.null(col))
      col <- sample(rainbow(max(y)))
    xlim <- c(1, length(n))
    z <-  x[g[i], n]
    if (leg) {
      ylab = "Expression"
    }
    else {
      ylab = NA
    }
    if (!is.null(name)) {
      main = name
    }
    else {
      main = NA
    }
    if (i == 1){
      if (is.null(ylim)) {
        plot(c(1, length(n)), c(min(z), max(z)), cex = 0, axes = FALSE,
             xlab = "", ylab = ylab, main = main, xlim = xlim)
      }
      else {
        plot(c(1, length(n)), c(min(z), max(z)), cex = 0, axes = T,
             xlab = "", ylab = ylab, main = main, xlim = xlim,
             ylim = ylim, xaxt='n')
      }}
    u <- 1:length(n)
    v <- as.vector(t(z))
    zc <- predict(loess(v ~ u, span = alpha))
    zc[zc < 0] <- 0.1
    lines(u, zc, lwd = lwd, col=set1[i])
    
  }
  width <- .45
  
  k <- 1:length(y)
  
  rect(xleft = k-width, ybottom = rep(-0.2, length(y)), xright =  k + width, ytop = rep(0, length(y)), col = col[y], border = NA)
  
  if ( samp == T) {
    if ( is.null(samppart)){
      stop("set samp part ")
    }
    rect(xleft = k-width, ybottom = rep(-0.5, length(y)), xright =  k + width, ytop = rep(-0.3, length(y)), col = sampcol[samppart], border = NA)
  }
  if (!leg)
    box(col = "white")
  else legend(legendpos, g, col=set1, pch=20)
}

plotsymbol_short <- function (d, types, subset = NULL, samples_col = NULL, cex = 0.5, 
                              fr = FALSE, um = FALSE, leg = TRUE, map = TRUE, cells=NULL) 
{
  subset2 <- subset
  if (is.null(subset)) 
    subset <- unique(types)
  h <- sort(unique(types)) %in% subset
  if (!is.null(subset)) {
    fp <- rep(FALSE, length(types))
    fp[types %in% subset] <- TRUE
    pchs <- rep(c(20,2), length(subset))
  }
  if (is.null(samples_col)) {
    samples_col <- rainbow(length(unique(types[fp])))
  }
  else {
    samples_col <- samples_col[h]
  }
  d <- d
  if (map) {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = cex, 
         pch = 20, col = "white")
    for (i in 1:length(unique(types[fp]))) {
      f <- types == sort(unique(types[fp]))[i]
      if (!is.null(cells)) {
        f2 <- rownames(d) %in% cells
        f <- f & f2 
      }
      if (!is.null(subset2)) {
        points(d[f, 1], d[f, 2], col = samples_col[i], pch = pchs[i], 
               cex = cex)
      }
      else{
        points(d[f, 1], d[f, 2], col = samples_col[i], pch = 20, 
               cex = cex)
      }
    }
  }
  else {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = 0, 
         pch = 20, col = "grey", xlim = c(min(d[, 1]), max(d[, 
                                                             1])), ylim = c(min(d[, 2]), max(d[, 2])))
  }
  if (leg) 
    if (!is.null(subset2)) {
      legend("bottomleft", legend = sort(unique(types[fp])), col = samples_col, 
             pch = pchs, cex = 0.75, bty = "n")
    }
  else{
    legend("bottomleft", legend = sort(unique(types[fp])), col = samples_col, 
           pch = 20, cex = 0.75, bty = "n")}
}

plotmap_short <- function (d, part, tp = 1, fr = FALSE, um = FALSE, 
                           cex = 0.5, leg=F, legendpos=NULL, fcol=NULL, medoids=NULL, clustsize=NULL) 
{ 
  d <- d
  part <- part
  if (!is.null(clustsize)) {
    uni_part <- sort(as.numeric(names(table(part)[table(part) >= clustsize])))
    index <- which(part %in% uni_part)
    d <- d[index,]
    part <- part[index]
    index_med <- sort(as.numeric(names(table(part)[table(part) >= clustsize])))
    medoids <- medoids[index_med]
  }
  
  else{
    uni_part <- sort(as.numeric(names(table(part))))
  }
  if (!is.null(fcol)) { fcol <- fcol}
  else{
    fcol <- rainbow(length(uni_part))}
  row.names(d) <- names(part)
  plot(d, xlab = "", ylab = "", cex = 0, axes = FALSE)
  for (i in 1:length(uni_part)) {
    if (sum(part == uni_part[i]) > 0) 
      points(d[part == uni_part[i], 1], d[part == uni_part[i], 2], col = adjustcolor(fcol[i], 
                                                                                     tp), pch = 20, cex = cex)
  }
  if (!is.null(medoids)) {
    for (i in 1:length(uni_part)) {
      if (sum(part == uni_part[i]) > 0) 
        points(d[medoids[i], 1], d[medoids[i], 
                                   2], col = adjustcolor(fcol[i], tp), pch = 20, 
               cex = 4)
      if (sum(part == uni_part[i]) > 0) 
        points(d[medoids[i], 1], d[medoids[i], 
                                   2], col = adjustcolor("white", tp), pch = 20, 
               cex = 3)
      if (sum(part == uni_part[i]) > 0) 
        text(d[medoids[i], 1], d[medoids[i], 
                                 2], uni_part[i], col = adjustcolor("black", tp), cex = 0.75, 
             font = 4)
    }
  }
  if(leg) {
    legend(legendpos, legend = sort(unique(part)), col = fcol, pch = 20, cex = 0.75, bty = "n" )
  }
}

get_cluster_composition_forest <- function(part, symbols=NULL, cluster_size, norm=T, color = NULL, map = T, leg=T, order_clus=NULL,ordertype=NULL, besides=T, enrichment=T, tosample=T){
  cpart <- part
  ### needs to be updated when input is partition of celltypes as factors
  if (class(part) %in% c("numeric", "factor")){
    cluster <- as.numeric(names(table(cpart)[table(cpart) >= cluster_size]))}
  else{
    cluster <- names(table(cpart)[table(cpart) >= cluster_size])
  }
  if (is.null(symbols)) {
    symbols <- sub("\\_.+", "", names(cpart))
  }
  types <- unique(symbols)
  types <- types[order(types)]
  if (!is.null(ordertype)) { types <- types[ordertype]}
  cpart <- cpart[cpart %in% cluster ]
  if (is.null(color)) {
    if ( tosample == T) {
      color <- rainbow(length(cluster))
    }
    color <- rainbow(length(types))
  }
  iels <- cbind(cluster=as.character(cpart),sample=symbols)
  rownames(iels) <- names(cpart)
  iels <- data.frame(iels)
  counts <- table(iels$sample,iels$cluster)
  if ( !is.null(ordertype)) { counts <- counts[ordertype,]}
  #
  if ( norm == T){
    counts <- counts/as.numeric(table(symbols))  
    if (tosample) {
      title_main <- paste0("Norm. ", "sample fraction of cluster")}
    else{
      title_main <- paste0("Norm. ", "cluster fraction of cluster")
    }
  }
  else{
    if (tosample){
      title_main <- "sample fraction of cluster"}
    else{
      title_main <- "cluster fraction of cluster"
    }
  }
  if ( !is.null(order_clus)) {
    rel_counts <- rel_counts[,as.character(order_clus)]
  }
  if (tosample == T) {
    rel_counts <- counts/apply(counts, 1, sum)
    if (class(part) %in% c("numeric", "factor")) {
      rel_counts <- rel_counts[,order(as.numeric(colnames(rel_counts)))]}
    rel_counts <- t(rel_counts)
    
  }
  else{
    rel_counts <- t(t(counts)/apply(counts,2,sum))
    if (class(part) %in% c("numeric", "factor")) {
      rel_counts <- rel_counts[,order(as.numeric(colnames(rel_counts)))]}
  }
  if (enrichment) {
    data_relcounts <- data.frame(rel_counts)
    if ( tosample == T) {
      colnames(data_relcounts)[1:2] <- c("Cluster", "Sample")
    }
    else{
      colnames(data_relcounts)[1:2] <- c("Sample", "Cluster")}
    
    
    library(tidyr)
    pop <- length(cpart)
    cluster <- sort(cluster)
    sample <- types
    pvals <- list()
    if ( tosample == T){
      for ( k in 1:length(cluster)) {
        clusteri <- cpart == cluster[k]
        clus_vec <- rep("nonclus", length(cpart))
        names(clus_vec) <- names(cpart)
        clus_vec[which(clusteri)] <- "tarclus"
        pvalues = c()
        pvalues2 <- c()
        for ( i in 1:length(sample)) {
          samperi <- grepl(sample[i], names(part))
          samp_vec <- rep("nonsamp", length(cpart))
          names(samp_vec) <- names(cpart)
          samp_vec[which(samperi)] <- "tarsamp"
          contigency <- table(clus_vec, samp_vec)
          fisher <- fisher.test(contigency, alternative = "greater")
          pv <- fisher$p.value
          pvalues <- append(pvalues, pv)        
          fisher2 <- fisher.test(contigency, alternative = "less")
          pv2 <- fisher2$p.value
          pvalues2 <- append(pvalues2, pv2)          
        }
        pval_data <- cbind(pvalues, pvalues2)
        rownames(pval_data) <- sample
        pvals[[as.character(cluster[k])]] <- pval_data
      }
    }
    else{
      for ( k in 1:length(sample)) {
        samperi <- grepl(sample[k], names(part))
        samp_vec <- rep("nonsamp", length(cpart))
        names(samp_vec) <- names(cpart)
        samp_vec[which(samperi)] <- "tarsamp"
        pvalues <- c()
        pvalues2 <- c()
        for ( i in 1:length(cluster)) {
          clusteri <- cpart == cluster[i]
          clus_vec <- rep("nonclus", length(cpart))
          names(clus_vec) <- names(cpart)
          clus_vec[which(clusteri)] <- "tarclus"
          contigency <- table(samp_vec, clus_vec)
          fisher <- fisher.test(contigency, alternative = "greater")
          pv <- fisher$p.value
          pvalues <- append(pvalues, pv)
          fisher2 <- fisher.test(contigency, alternative = "less")
          pv2 <- fisher2$p.value
          pvalues2 <- append(pvalues2, pv2)
        }
        pval_data <- cbind(pvalues, pvalues2)
        rownames(pval_data) <- as.character(cluster)
        pvals[[sample[k]]] <- pval_data
      }}
    pvals <- lapply(pvals, function(x){
      apply(x, 2, function(y){
        p.adjust(y, method="bonferroni")
      })
    })
    
    pvals2 <- lapply(pvals, function(x, tosample) {
      x <- data.frame(x)
      if ( tosample==T){ x <- cbind(x, Sample=rownames(x))}
      else{ x <- cbind(x, Cluster=rownames(x))}
      x
    }, tosample=tosample)
    
    pvals3 <- lapply(seq_along(pvals2), function(x, pvlist, tosample) {
      y <- pvals2[[x]]
      nam <- names(pvlist)[x]
      if (tosample==T){ y <- cbind(y, Cluster=rep(nam, nrow(y)))}
      else{y <- cbind(y, Sample=rep(nam, nrow(y)))}
      return(y)
    }, pvlist=pvals2, tosample=tosample)
    
    pvals4 <- do.call(rbind, pvals3)
    if (tosample==T) {pvals4 <- pvals4[order(pvals4$Sample),]}
    else{
      if (class(part) %in% c("numeric", "factor")) {
        pvals4 <- pvals4[order(as.numeric(pvals4$Cluster)),]}
      else{
        pvals4 <- pvals4[order(pvals4$Cluster),]
      }
    }
    
    data_relcounts$Cluster <- as.character(data_relcounts$Cluster)
    data_relcounts$Sample <- as.character(data_relcounts$Sample)
    
    merged_data <- merge(data_relcounts, pvals4, by=c("Cluster", "Sample"))
    merged_data <- cbind(merged_data, differential=rep("n.s.", nrow(merged_data)))
    merged_data[which(merged_data$pvalues < 0.05),"differential"] <- "hi"
    merged_data[which(merged_data$pvalues2 < 0.05),"differential"] <- "lo"
    if (class(part) == "numeric"){
      merged_data$Cluster <- factor(merged_data$Cluster, levels = c(unique(merged_data$Cluster)[order(as.numeric(unique(merged_data$Cluster)))]))}
    if ( tosample==T){
      if (norm) {
        p2 <- ggplot(merged_data, aes(fill=Cluster, y=Freq, x=Sample, label= differential)) + 
          geom_bar(stat="identity") + ggtitle(expression(atop("Norm. cluster fraction of sample", atop("high/low/n.s. = bonferroni adjusted p.value <0.05 of enrichment")))) +
          scale_fill_manual(values=color) + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
      else{
        p2 <- ggplot(merged_data, aes(fill=Cluster, y=Freq, x=Sample, label= differential)) + 
          geom_bar(stat="identity") + ggtitle(expression(atop("cluster fraction of sample", atop("high/low/n.s. = bonferroni adjusted p.value <0.05 of enrichment")))) +
          scale_fill_manual(values=color) + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      }
    }
    else{
      if (norm){
        p2 <- ggplot(merged_data, aes(fill=Sample, y=Freq, x=Cluster, label= differential)) + 
          geom_bar(stat="identity") + ggtitle(expression(atop("Norm. sample fraction of cluster", atop("high/low/n.s. = bonferroni adjusted p.value <0.05 of enrichment")))) +
          scale_fill_manual(values=color) + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))}
      else{
        p2 <- ggplot(merged_data, aes(fill=Sample, y=Freq, x=Cluster, label= differential)) + 
          geom_bar(stat="identity") + ggtitle(expression(atop("sample fraction of cluster", atop("high/low/n.s. = bonferroni adjusted p.value <0.05 of enrichment")))) +
          scale_fill_manual(values=color) + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      }
      
    }
    print(p2)
    
  }
  else{
    if (besides==T) {
      if(map == F) {
        barplot(rel_counts, beside=T, col = NA, border = NA, axes = FALSE, axisnames = F)
      }
      else{
        barplot(rel_counts, beside=T,border = NA, main=title_main, ylab="% of cluster",
                xlab="", col=color, names.arg=as.character(colnames(rel_counts)), cex.names=1, las=2)
      }
      if (leg==T) {
        legend( "topright", pch=20, bty="n", cex=1, legend=types, col=color) 
      }
    }
    else{
      if (map == F) {
        barplot(rel_counts, col = NA, border = NA, axes = FALSE, axisnames = F)
      }
      else{
        barplot(rel_counts, main=title_main, ylab="% of cluster",
                xlab="", col=color, names.arg=as.character(colnames(rel_counts)), cex.names=1, las=2)
      }
      if (leg==T) {
        legend( "topright", pch=20, bty="n", cex=1, legend=types, col=color) 
      }
    }
  }
}

get_agg_sd <- function(sd, cent=F, batch=NULL) {
  if (cent==F) {
    sd_filt <- sd[,-((ncol(sd)-2):ncol(sd))]
    sample <- sd[,"batch"]
  }
  
  else{
    sd_filt <- sd
    sample = batch
  }
  comb_sd_num <- apply(sd_filt, 2, as.numeric)
  rownames(comb_sd_num) <- rownames(sd_filt)
  
  agg3 <- aggregate(comb_sd_num, list(sample), mean)
  
  rownames(agg3) <- agg3[,1]
  agg3 <- agg3[,-1]
  agg3
}
plot_sd_comp <- function(sd1, sd2, nam1, nam2,fcol, legendpos="topleft", path) {
  if(identical(rownames(sd1), rownames(sd2))) {
    nam <- paste0(nam1, "_vs_", nam2, ".pdf")
    pdf(file.path(path, nam))
    plot(sd1[,"RPE"], sd2[,"RPE"], xlab=nam1,ylab=nam2,pch=20, cex=1.5, col=fcol)
    legend(legendpos, rownames(sd1), fill=fcol, bty="n")
    dev.off()
  }
  else{
    cat("not the same samples of sd1 versus sd2")
  }
}

