
## function to plot heatmap clustering of genes and samples from seurat object
#' DotPlot(object = pbmc_small, features = cd_genes)
#'
plotExp <- function(
  object,
  features,
  idents = NULL,
  col.min = -3,
  col.max = 2,
  scale = TRUE,
 group.by = NULL
) {

## extract data to plot by cells 
cells <- unlist(x = SeuratObject::CellsByIdentities(object = object, idents = idents))
data.features <- FetchData(object = object, vars = features, cells = cells)
  
 ## add group id as a column

  data.features$id <- if (is.null(x = group.by)) {
    Idents(object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)

## caclulate the avg.exp & percentage of cells expressed
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  
names(x = data.plot) <- unique(x = data.features$id)

## structure data to plot with ggplot2

# convert into data.frame
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  
 data.plot <- do.call(what = 'rbind', args = data.plot)

#scale the avg exp
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  
 data.plot$avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
 
 data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )

 data.plot$id <- factor(data.plot$id, levels=rev(cells.type.order))
 data.plot <- data.plot %>% filter(pct.exp > 0, avg.exp.scaled > 0 )
 
 plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', fill = 'avg.exp.scaled'), colour="black",pch=21) +
#    scale_fill_distiller(palette = "RdBu")+
    scale_fill_gradient(low="white", high="darkblue", name="Mean expression \n in group",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = 15, barheight = 1))+
    theme_bw()+
    theme(text=element_text(size=22, color="black"),
                        axis.text.x =element_text(size=20, color="black", angle = 90, hjust = 1, vjust = 1),
                        axis.text.y =element_text(size=18, color="black"),
                         legend.position="top",
                        panel.grid.major = element_line(size = 0.75),
                        panel.grid.minor = element_line(size = 0.75),
                        panel.border = element_rect(colour = "black", fill=NA, size=1.5))+
    labs(x="", y="")+
    scale_size(range=c(1,7))+
#    scale_radius()+
    guides(size = guide_legend(title = 'Fraction of cells \n in group'))


return(plot)
	
}



# helpping functions
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
