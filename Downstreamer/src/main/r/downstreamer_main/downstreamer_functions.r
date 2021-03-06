library(ggplot2)
library(gridExtra)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(data.table)

# Function definitions
# ------------------------------------------------------ 
# NPG friendly ggplot theme
theme.nature <- function(p, base_size = 11, base_family = "ArialMT") {
  p <- p + theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black", size=0.5),
          axis.ticks = element_line(size=0.5),
          axis.text = element_text(size=base_size, family="ArialMT", face="plain"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=base_size, family = "ArialMT", face="plain"),
          complete = TRUE,
          plot.title = element_text(hjust=0.5))
  return(p)
}

# ------------------------------------------------------ 
# Plot depict2 scatterplot by enrichment score
make.tsne.plot <- function(data, trait, x="Annotation1", y="Annotation2", colour="Enrichment.Z.score", limits=NA) {
  data   <- data[order(abs(data[,colour])),]
  
  if (is.na(limits)) {
    limits <- c(-max(abs(data[,colour])), max(abs(data[,colour])))
  } else {
    data[data[,colour] < limits[1]+0.0001 ,colour] <- limits[1]
    data[data[,colour] > limits[2]-0.0001 ,colour] <- limits[2]
    
  }
  
  # low  <- colorRampPalette(c('firebrick3', '#F6DCDC'))
  # mid  <- colorRampPalette(c("#F6DCDC", "#DDE6EE"))
  # high <- colorRampPalette(c('#DDE6EE', 'dodgerblue4'))
  # bins= 100
  # breaks = seq(-15, 15, length.out = bins)
  # colRamp <- c(low((bins-10)/2), mid(10), high((bins-10)/2))
  # cols <- colRamp[as.numeric(cut(data[,colour], breaks = breaks))]
  # 
  # pdf(width=7.5, height=7.5, file="~/Desktop/")
  # 
  # 
  # layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
  # plot(data$Annotation1, data$Annotation2, pch = 20, col=adjustcolor(cols, alpha.f = 0.8), bty='n', xlab="t-SNE component 1", ylab="t-SNE compenent 2")
  # 
  # legend_image <- as.raster(matrix(colRamp, ncol=1))
  # plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Z-score')
  # text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
  # rasterImage(legend_image, 0, 0, 1,1)
  # 
  #  legend("topright",legend=c("15", "0", "-15"),fill=c("dodgerblue4", "grey", "firebrick3"), title="Z-score", border=F, bty="n")
  # #legend.scale(c(-15, 15), col=cols)
  # 
   p <- ggplot(aes(x=as.numeric(data[,x]),
                   y=as.numeric(data[,y]),
                   colour=data[,colour]),
               data=data) +
     geom_point() +
     labs(title=trait) +
     xlab(x) +
     ylab(y) +
     scale_colour_gradient2(name="Z-score", limits=limits, low="firebrick3", mid="white", high="dodgerblue4")
  
  #p <- ggplot(aes(x=as.numeric(data[,x]), y=as.numeric(data[,y])), data=data) +
  #  labs(title=trait) + xlab(x) + ylab(y)
  
  #p <- ggplot(aes(x=as.numeric(data[,x]), y=as.numeric(data[,y]), colour=data[,colour]), data=data) +
  #  geom_point() +
  #  labs(title=trait) + xlab(x) + ylab(y) +
  #  scale_colour_gradient(low=adjustcolor("red", alpha.f = 0.5),
  #                        high=adjustcolor("blue", alpha.f = 0.5),
  #                        limits=c(-15, 15),
  #                        name="Z-score")
  
  p <- theme.nature(p)
  return(p)
}

# ------------------------------------------------------ 
read.depict2 <- function(path, potential_traits=NULL) {
  if (is.null(potential_traits)) {
    potential_traits <- excel_sheets(path)
    potential_traits <- potential_traits[grep("Overview", potential_traits, invert=T)]
  }

  output <- list()
  for (sheet in potential_traits) {
    tmp <- tryCatch({data.frame(read_excel(path, sheet=sheet, col_types ="guess", trim_ws = T), stringsAsFactors=F)},
                    error=function(a){return(NA)},
                    warn=function(a){return(NA)})
    
    if(!is.na(tmp)) {
      
      for (i in 1:ncol(tmp)) {
        if (class(tmp[,i]) == "character"){
          tmp[,i] <- type.convert(tmp[,i], as.is=T)
          
        }
      }
      
      rownames(tmp) <- tmp[,1]
      output[[sheet]] <- tmp
    }
  }
  
  return(output)
}

# ------------------------------------------------------ 
# Read a batch of downstreamer results
read.downstreamer.batch <- function(main.downstreamer.output.path, potential_traits=NULL, USE.CACHE=F) {

  if (!dir.exists("data")) {
    dir.create("data")
  }

  # Loads enrichment results
  if (USE.CACHE & file.exists("data/downstreamer_results_cache.RData")) {
    cat("[INFO] Loading downstreamer results from cache\n")
    load("data/downstreamer_results_cache.RData")
  } else {
    cat("[INFO] Loading downstreamer results from files\n")
    files    <- list.files(main.downstreamer.output.path, pattern="*\\_enrichtments\\_exHla\\.xlsx", full.names = T)
    datasets <- list()
    for (file in files) {
      name <- gsub("\\_hg19\\_enrichtments\\_exHla\\.xlsx", "", basename(file))
      name <- gsub("\\_hg19\\_enrichtments\\_exHla\\_1\\.xlsx", "", name)
      name <- gsub("\\_enrichtments\\_exHla\\.xlsx", "", name)
      name <- gsub("\\_hg19\\.txt\\_exHla\\.xlsx", "", name)

      datasets[[name]] <- read.depict2(file, potential_traits = potential_traits)
    }
    cat("[INFO] Done, saving cache for future use\n")
    save(datasets, file="data/downstreamer_results_cache.RData")
  }

  return(datasets)
}

# ------------------------------------------------------
# Convert a list of downstreamer enrichments into a matrix
make.zscore.matrix <- function(datasets, trait="Coregulation", collumn="Enrichment.Z.score") {
  out <- matrix()
  i=0
  for (dataset in names(datasets)) {
    tmp <- datasets[[dataset]][[trait]]
    if (i==0) {
      out <- matrix(tmp[,collumn])
      rownames(out) <- tmp[,1]
    } else {
      out <- cbind(out, tmp[rownames(out), collumn])
    }
    i <- i+1
  }
  colnames(out) <- names(datasets)
  
  return(out)
}

# ------------------------------------------------------ 
read.enrichments <- function(files, column=1, trim.colnames=NULL) {
  out <- matrix()
  i=0
  for (dataset in files) {
    tmp <- read.table(dataset, stringsAsFactors = F, row.names = 1, header=T, sep="\t")
    if (i==0) {
      out <- as.matrix(tmp[,column, drop=F])
    } else {
      out <- cbind(out, tmp[rownames(out), column, drop=F])
    }
    i <- i+1
  }
  colnames(out) <- make.names(basename(files), unique=T)
  
  if (!is.null(trim.colnames)) {
    colnames(out) <- gsub(trim.colnames, "", colnames(out))
  }
  
  return(out)
}

# ------------------------------------------------------ 
read.genep.excel <- function(files, column=6) {
  out <- matrix()
  i=0
  
  rnames <- c()
  cnames <- c()
  
  for (dataset in files) {
    tmp <- data.frame(read_excel(dataset, sheet=1, col_types ="guess", trim_ws = T), stringsAsFactors = F)
    rownames(tmp) <- tmp[,1]
    if (i==0) {
      out <- as.matrix(tmp[,column:ncol(tmp), drop=F])
      rownames(out) <- tmp[,1]
      rnames <- tmp[,1]
    } else {
      out <- cbind(out, tmp[rownames(out), column:ncol(tmp), drop=F])
    }
    
    if (ncol(tmp) - column > 0) {
      cnames <- c(cnames, colnames(tmp)[column:ncol(tmp)])
    } else {
      cnames <- c(cnames, basename(dataset))
    }
    
    i <- i+1
  }

  out <- apply(out, 2, as.numeric)
  colnames(out) <- cnames
  rownames(out) <- rnames
  
  return(out)
}

# ------------------------------------------------------ 
read.enrichments.fread <- function(files, column=1) {
  out <- matrix()
  i=0
  for (dataset in files) {
    tmp <- fread(dataset, stringsAsFactors = F, header=T, sep="\t", data.table=F)
    rownames(tmp) <- make.names(tmp[,1], unique=T)
    tmp <- tmp[,-1]
    if (i==0) {
      out <- as.matrix(tmp[,column, drop=F])
    } else {
      out <- cbind(out, tmp[rownames(out), column, drop=F])
    }
    i <- i+1
  }
  colnames(out) <- basename(files)
  
  return(out)
}

# ------------------------------------------------------ 
cor.test.p <- function(x, method){
  FUN <- function(x, y) cor.test(x, y, method=method, use="complete.obs")[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

# ------------------------------------------------------ 
make.correlation.heatmap <- function(matrix, trait) {
  c <- cor(matrix, method="spearman", use="pairwise.complete.obs")
  #p <- cor.test.p(matrix, method="spearman")
  
  #c[p > 0.05] <- 0
  
  hm(c, main=trait)
}

# ------------------------------------------------------ 
hm <- function(data, cellwidth=12, cellheight=12, limit=NULL, ...) {
  break.list <- seq(-max(abs(data)), max(abs(data)), by=max(abs(data))/100)
  pheatmap(data,
           breaks=break.list,
           col=colorRampPalette(rev(brewer.pal(n=7, name ="RdBu")))(length(break.list)),
           cellwidth=cellwidth,
           cellheight=cellheight,
           ...)
}

# ------------------------------------------------------ 
simple.qq.plot <- function (observedPValues) {
  plot(-log10(1:length(observedPValues)/length(observedPValues)), 
       -log10(sort(observedPValues)))
  abline(0, 1, col = "red")
}

#----------------------------------------------------------------------------------------
fancy.qq.plot <- function (y, main="", col="black", highlight.col="blue") {
  y <- na.omit(y)
  x <- -log10(1:length(y)/length(y))
  y <- -log10(sort(y))
  
  thresh <- -log10(0.05/length(y))
  upper <- max(c(x, y)) + 1
  
  p <- ggplot(data.frame(x=x, y=y, col=y > thresh), mapping=aes(x=x, y=y, col=col)) +
    geom_point() +
    geom_abline(slope=1, intercept=0, col="black") +
    geom_hline(yintercept=thresh, lty=2, col="grey") +
    xlab("-log10(expected p-value)") +
    ylab("-log10(observed p-value)") +
    ggtitle(main) +
    ylim(c(0, upper)) +
    xlim(c(0, upper)) +
    coord_fixed() +
    annotate("text", y=01, x=upper-2, label=paste0("N=", length(y)), col="grey") +
    scale_color_manual(values=c(`TRUE`=highlight.col, `FALSE`=col), name="Bonferroni significant")
  
  return(p)
}

#----------------------------------------------------------------------------------------
xy.plot.pvalue.colored <- function(auc.1, auc.pval.1, auc.2, auc.pval.2, xlab="X", ylab="Y", main=NULL) {
  auc.pval.1[auc.1 == 0] <- 1
  auc.pval.2[auc.2 == 0] <- 1
  
  df.plot <- data.frame(auc.1=auc.1,
                        auc.2=auc.2,
                        signif.1 = auc.pval.1 < (0.05 / length(auc.1)), 
                        signif.2 = auc.pval.2 < (0.05 / length(auc.2)))
  
  df.plot$signif         <- sum(df.plot$signif.1 + df.plot$signif.2)
  df.plot$signif.both    <- (df.plot$signif.1 + df.plot$signif.2) == 2
  df.plot$signif.either  <- (df.plot$signif.1 + df.plot$signif.2) > 0
  df.plot <- df.plot[order(df.plot$signif.either),]
  
  lims <- c(min(c(auc.1, auc.2), na.rm=T), max(c(auc.1, auc.2), na.rm=T))  

  p <- ggplot(data=df.plot, mapping=aes(x=auc.1, y=auc.2)) +
    geom_point(alpha=0.75, mapping=aes(col=signif.either)) +
    geom_abline(slope=1, intercept=0, col="grey", lty=2) +
    coord_fixed() +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(main) + 
    scale_color_manual(values=c(`FALSE`="#2c6c70", `TRUE`="#0ae4f2")) +
    geom_smooth(method="lm") +
    xlim(lims) +
    ylim(lims)
  
  return(theme.nature(p))
}



