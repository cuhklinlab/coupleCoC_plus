# plot UMAP for the raw data
plot.umap = function(x, labels, main, colors, pos = "topright", pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="", cex.main=1, cex.legend=0.85) {
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  }

  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F, xlab = "UMAP1")
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)


  labels.u = unique(labels)
  legend.pos = pos
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }

  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

# functions for plotting the heatmaps
clu_sep <- function(CZ_best){
  v <- c()
  for (i in 1:length(unique(CZ_best))){
    v = c(v,which(CZ_best==i))
  }
  return (v)
}

line_sep <- function(CZ_best){
  v <- c()
  v[1] <- length(which(CZ_best==1))
  for (i in 2:length(unique(CZ_best))){
    v[i] = v[i-1]+length(which(CZ_best==i))
  }
  return (v)
}

clu_num <- function(CZ_best){
  v <- c()
  for (i in 1:length(unique(CZ_best))){
    v = c(v,length(which(CZ_best==i)))
  }
  return (v)
}

# strenghten the signals
strong_signal <- function(raw_data,CX_best,rowsep,n_pick){
  #raw_data = X_raw[cell_order_X,gene_order_Z];rowsep = rowsep_x; n_pick = 15;
  #raw_data = YY;rowsep = rowsep_y; n_pick = 15;CX_best=CY_best;
  n_cluster = length(unique(CX_best))
  b_vec = c(0,cumsum(rowsep))
  for (i in 1:n_cluster){
    # for the j-th block
    j = (b_vec[i]+1):b_vec[i+1]
    data = raw_data[j,]
    for(k in 1:dim(data)[1]){
      # for the k-th sample in the j-th block
      select = sample(1:dim(data)[1], size=n_pick, replace=F)
      data[k,] = colMeans(data[select,],na.rm = TRUE)#colSums(data[select,])/n_pick
    }
    raw_data[j,]=data
  }
  return(raw_data)
}

#reorder the feature sequence based on the vertical signal strength
feature_signal_reorder<-function(colsep,X){
  n_cluster = length(colsep)
  z_vec = c(0,colsep)
  for(i in 1:n_cluster){
    # for the j-th block
    j = (z_vec[i]+1):z_vec[i+1]
    data = X[,j]
    jj = order(as.vector(colSums(data,na.rm = TRUE)))
    X[,j] = data[,jj]
  }
  return(X)
}


heatmap_fun <- function(X, scaleyellowred, colsep, rowsep_x){
  heatmap.2(X, dendrogram="none", col = scaleyellowred, na.color="yellow",
            labRow = "", labCol = "",
            ## row and column order
            Rowv=FALSE, Colv = FALSE,
            ## row and column side color
            #ColSideColors = colcolssubS
            #RowSideColors = cols,
            margins = c(2, 2),  trace = "none",
            ## color key
            key=T, key.xlab="log2(normalized read count + 1)", key.title="", #keysize=1.2,
            density.info = "none",
            #main="threshold by 3"
            # plot labels
            #main = "K562 and HL60 scATAC_seq data after clustering by coupleCoC",
            xlab = "features",
            ylab = "cells",
            sepcolor="blue",
            lwd = 3,
            colsep=colsep,
            rowsep=rowsep_x,
  )
}


heatmap_fun_true <- function(X, scaleyellowred){
  heatmap.2(X, dendrogram="none", col = scaleyellowred, na.color="yellow",
            labRow = "", labCol = "",
            ## row and column order
            Rowv=FALSE, Colv = FALSE,
            ## row and column side color
            #ColSideColors = colcolssubS
            #RowSideColors = cols,
            margins = c(2, 2),  trace = "none",
            ## color key
            key=T, key.xlab="log2(normalized read count + 1)", key.title="", #keysize=1.2,
            density.info = "none",
            #main="threshold by 3"
            # plot labels
            #main = "K562 and HL60 scATAC_seq data after clustering by coupleCoC",
            xlab = "features",
            ylab = "cells",
            sepcolor="purple",
  )
}

#re-order the label to ensure

reorder_label <- function(CX_best,old_label){
  CX0 <- CX_best
  #old_label<-c(2,1)
  labelen <-length(unique(old_label))
  new_label<-1:labelen
  for (i in 1:labelen){
    CX0[which(CX0==old_label[i])] = labelen+i
  }
  new_label = 1
  for (j in (labelen+1):(2*labelen)){
    CX0[which(CX0==j)] = new_label
    new_label = new_label + 1
  }
  return(CX0)
}
