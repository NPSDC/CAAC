replace.symbol <- function(sym.list, replace.with, to.replace)
{
  sapply(sym.list, function(x)
    {
    temp <- unlist(strsplit(x, to.replace, fixed = T))
    paste(temp, collapse = '-')
  })
}

get.imp.genes <- function(method, assay, number)
{
  if(method == 1) ###most expression
    return(order(rowMeans(assay), decreasing = T)[1:number])
  else if(method == 2)
  {
    mads = apply(assay, 1, mad)
    return(order(mads, decreasing = T)[1:number])
  }
  else if(method == 3)
  {
    library(matrixStats)
    rv <- rowVars(assay)
    return(order(rv, decreasing=TRUE)[1:number])
  }
    
}

create.list.imp.genes <- function(df.list, method, number)
{
  top.genes = list()
  top.genes[['rld']] = get.imp.genes(method, assay(df.list[['rld']]), number)
  top.genes[['nt']] = get.imp.genes(1, assay(df.list[['nt']]), number)
  top.genes[['vs']] = get.imp.genes(1, assay(df.list[['vs']]), number)
  top.genes[['fpqm']] = get.imp.genes(1, df.list[['fpqm']], number)
  top.genes[['fpqm_log']] = get.imp.genes(1, df.list[['fpqm_log']], number)
  return(top.genes)
  
}

##Tweeked original multiplot and plotPCA functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  if(length(...) == 1)
  {
    print('aye')
    plots <- c(list(...), plotlist)
  }
  else
    plots <- c(..., plotlist)
  #print(plots)
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE, title, sample.info = NULL)
{
  library(ggplot2)
  library(matrixStats)
  # calculate the variance for each gene
  if(typeof(object) == 'S4')
  {
    #rv <- rowVars(assay(object))
  
    # select the ntop genes by variance
    #select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
    
    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
      factor(apply( intgroup.df, 1, paste, collapse=" : "))
    } else {
      colData(object)[[intgroup]]
    }
  }
  else
  {
    
    # perform a PCA on the data in assay(x) for the selected genes
    
    object = na.omit(object)
    pca <- prcomp(t(object))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    if (!all(intgroup %in% names(sample.info))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(sample.info[, intgroup, drop=FALSE])
    
    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
      factor(apply( intgroup.df, 1, paste, collapse=" : "))
    } else {
      sample.info[[intgroup]]
    }
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    ggtitle(title)+
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() 
}

get.list.plotPCA <- function(arg, assay, gene.list, main_string, type = 1, sample.info = NULL, intgroup = 'type')
{
  if(length(arg) == 0& length(gene.list) == 0)
    li <-plotPCA(assay, intgroup = 'type', title = main_string)
  else
  {
    if(type == 1)
    {
      li <- lapply(arg, function(x)
        {
        a = plotPCA(assay[gene.list[[x]],], intgroup = intgroup, 
                title = paste(main_string,as.character(x), sep = '_'), sample.info = sample.info)})
      
        print(a)
    }
    else if(type == 2)
    {
      li <- lapply(arg, function(x)
      {
        plotPCA(assay[gene.list[1:x],], intgroup = intgroup, 
                title = paste(main_string,as.character(x), sep = '_'), sample.info = sample.info)})
    }
      }
  return(li)
}