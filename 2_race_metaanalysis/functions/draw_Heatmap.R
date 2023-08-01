# adapted from the bioconductor package "DExMA"
# added more parameters to customize heatmap further

draw_Heatmap <- function(objectMA, resMA,
                        typeMethod = c("FEM", "REM", "maxP", "minP","Fisher", "Stouffer", "ACAT"),
                        scaling = c("zscor","rscale","swr","none"),
                        regulation = c("all", "up","down"),
                        breaks = c(-2,2),
                        fdrSig = 0.05,
                        numSig,
                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                        na_col = "black",
                        legend = TRUE,
                        control,
                        case,
                        title)
{
  typeMethod <- match.arg(typeMethod)
  scaling <- match.arg(scaling)
  regulation <- match.arg(regulation)
  heatmap.var = "none"

  if(!is.null(breaks)){
    bks = seq(breaks[1],breaks[2],length.out = 100)
  } else{
    bks = NULL
  }

  ## Get features based on Meta-analysis results
  ## Case of Effects size
  if(typeMethod == "FEM" | typeMethod == "REM"){
    sig.genes <- .effectSig(resMA, regulation, numSig, fdrSig)
  }

  ## Case of P-value combination
  if(typeMethod == "Fisher" |typeMethod == "Stouffer" | typeMethod == "maxP" |
     typeMethod == "minP" | typeMethod == "ACAT"){
    sig.genes <- .pvalueSig(resMA, regulation, numSig, fdrSig)
  }

  ## Create a matrix with all datasets
  all.cl <- NULL
  all.colnames <- NULL
  sizes <- NULL
  for(i in seq_len(length(objectMA))){
    sizes <- c(sizes,rep(names(objectMA[i]), length(objectMA[[i]][[2]])))
    all.cl <- c(all.cl,objectMA[[i]][[2]])
    colnames(objectMA[[i]][[1]]) <- paste0(colnames(objectMA[[i]][[1]]),
                                           ".", names(objectMA[i]), 
                                           sep="")
    all.colnames <- c(all.colnames,colnames(objectMA[[i]][[1]]))
  }
  exp.ALL <- matrix(data=NA, ncol=length(all.colnames),nrow=length(sig.genes))
  rownames(exp.ALL) <- sig.genes
  colnames(exp.ALL) <- all.colnames

  ## Scaling datasets
  switch(scaling,
         rscale={exp.ALL<- .rscale(objectMA, sig.genes, exp.ALL)
         bks=NULL},
         zscor={exp.ALL<- .zscor(objectMA, sig.genes, exp.ALL)},
         none={heatmap.var="row"
         exp.ALL<- .none(objectMA, sig.genes, exp.ALL)},
         swr={heatmap.var="row"
         exp.ALL<- .swr(objectMA, resMA, sig.genes, exp.ALL)}
  )

  ## Heatmap
  ann <- data.frame(variable = ifelse(all.cl == 0, control, case),
                    dataset=sizes, row.names=colnames(exp.ALL))
  if(numSig>60){show.rownames=FALSE}
  else{show.rownames=TRUE}
  pheatmap(exp.ALL[,order(all.cl, decreasing=TRUE)], annotation_col = ann,
           cluster_cols=FALSE, cluster_rows=FALSE, scale=heatmap.var, breaks=bks,
           show_colnames=FALSE,show_rownames=show.rownames, na_col = na_col, 
           color = color, legend = legend, main = title)
  return(exp.ALL)
}

################################################################################

# Significant genes in Effects sizes results
.effectSig <- function(resMA, regulation, numSig, fdrSig){
  if(regulation=="up"){resMA <- resMA[rownames(resMA)[resMA$Zval>0],]}
  if(regulation=="down"){resMA <- resMA[rownames(resMA)[resMA$Zval<0],]}
  sig.genes <- rownames(resMA)[resMA$FDR<=fdrSig]
  ord <- resMA[sig.genes,]
  sig.genes <- sig.genes[order(ord$FDR,decreasing = FALSE)]
  if(length(sig.genes) < numSig){numSig <- length(sig.genes)}
  sig.genes <- sig.genes[seq_len(numSig)]
  ord <- resMA[sig.genes,]
  sig.genes <- sig.genes[order(ord$Zval,decreasing = TRUE)]
  return(sig.genes)
}

#Significant genes in p-values methods results
.pvalueSig <- function(resMA, regulation, numSig, fdrSig){
  if(regulation=="up"){resMA <- resMA[rownames(resMA)[resMA$AveFC>0],]}
  if(regulation=="down"){resMA <- resMA[rownames(resMA)[resMA$AveFC<0],]}
  sig.genes <- rownames(resMA)[resMA$FDR<=fdrSig]
  ord <- resMA[sig.genes,]
  sig.genes <- sig.genes[order(ord$FDR,decreasing = FALSE)]
  if(length(sig.genes) < numSig){numSig <- length(sig.genes)}
  sig.genes <- sig.genes[seq_len(numSig)]
  ord <- resMA[sig.genes,]
  sig.genes <- sig.genes[order(ord$AveFC,decreasing = TRUE)]
}

##############################################################################

# rscale function
.rscale <- function(objectMA, sig.genes, exp.ALL){
  print("scaling using rescale function...")
  for(set in seq_len(length(objectMA))){
    temp<-objectMA[[set]][[1]][intersect(sig.genes,
                                         rownames(
                                           objectMA[[set]][[1]]
                                         )),]
    
    for(gene in seq_len(nrow(temp))){
      temp[gene,] <- rescale(x=temp[gene,],from=c(min(temp[gene,]), 
                                                  max(temp[gene,])),
                             to=c(-1,1))
    }
    exp.ALL[rownames(temp),colnames(temp)] <- temp
  }
  return(exp.ALL)
}

# scaling zscor function
.zscor <- function(objectMA, sig.genes, exp.ALL){
  print("scaling using z-score...")
  for(set in seq_len(length(objectMA))){
    temp<-objectMA[[set]][[1]][intersect(
      sig.genes,
      rownames(objectMA[[set]][[1]])),]
    
    for(gene in seq_len(nrow(temp))){
      temp[gene,] <- (temp[gene,]-mean(temp[gene,],
                                       na.rm=TRUE))/sd(temp[gene,], 
                                                       na.rm = TRUE)
    }
    exp.ALL[rownames(temp),colnames(temp)] <- temp
  }
  return(exp.ALL)
}

# scaling with reference scaling function
.swr <- function(objectMA,resMA,sig.genes, exp.ALL){
  allgenes <- rownames(resMA)
  for (study in seq_len(length(objectMA))){
    studygenes <- rownames(objectMA[[study]][[1]])
    union <- unique(allgenes %in% studygenes)
    if(all(union)){
      ref <- names(objectMA[study])
      break}
    else{ref <- NULL}
  }
  if(is.null(ref)){
    stop("scaling with reference method can not be used. Please use 
        another method")
  }
  else{
    print("scaling with reference...")
    refH<-objectMA[[ref]][[1]][intersect(sig.genes,
                                         rownames(objectMA[[ref]][[1]])),
                               objectMA[[ref]][[2]]==0]
    refC<-objectMA[[ref]][[1]][intersect(sig.genes,
                                         rownames(objectMA[[ref]][[1]])),
                               objectMA[[ref]][[2]]==1]
    exp.ALL[rownames(refH),colnames(refH)] <- refH
    exp.ALL[rownames(refC),colnames(refC)] <- refC
    for(set in names(objectMA)){
      if(set!=ref){
        tempH<-objectMA[[set]][[1]][intersect(
          sig.genes,
          rownames(objectMA[[set]][[1]])),
          objectMA[[set]][[2]]==0]
        tempC<-objectMA[[set]][[1]][intersect(
          sig.genes, rownames(objectMA[[set]][[1]])),
          objectMA[[set]][[2]]==1]
        for(gene in seq_len(nrow(tempH))){
          tempH[gene,]<-(tempH[gene,]*(sd(refH[gene,],
                                          na.rm=TRUE)/sd(tempH[gene,],
                                                         na.rm=TRUE)))-
            (((mean(tempH[gene,],na.rm=TRUE)*
                 (sd(refH[gene,],na.rm=TRUE)/sd(tempH[gene,],
                                                na.rm=TRUE)))-
                mean(refH[gene,],na.rm=TRUE)))
          tempC[gene,]<-(tempC[gene,]*
                           (sd(refC[gene,],na.rm=TRUE)/
                              sd(tempC[gene,],na.rm=TRUE)))-
            (((mean(tempC[gene,],na.rm=TRUE)*
                 (sd(refC[gene,],na.rm=TRUE)/sd(tempC[gene,],
                                                na.rm=TRUE)))-
                mean(refC[gene,],na.rm=TRUE)))
        }
        exp.ALL[rownames(tempH),colnames(tempH)] <- tempH
        exp.ALL[rownames(tempC),colnames(tempC)] <- tempC
      }
    }
  }
  return(exp.ALL)
}

# none scaling function
.none <- function(objectMA, sig.genes, exp.ALL){
  for(set in seq_len(length(objectMA))){
    temp<-objectMA[[set]][[1]][intersect(
      sig.genes,rownames(objectMA[[set]][[1]])),]
    exp.ALL[rownames(temp),colnames(temp)] <- temp
  }
  return(exp.ALL)
}