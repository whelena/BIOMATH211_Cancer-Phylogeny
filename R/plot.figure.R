plot.figure <- function(
    tree,
    DT,
    fname = NULL,
    ...
    ) {

    DF <- data.frame(DT[, grep('snv', names(DT), value = TRUE), with = FALSE]);
    rownames(DF) <- DT$cell.id;
    phy <- ggtree(as.phylo(tree)) + geom_tippoint(aes(colour = true.col), color = tree$true.col, size = 3);
    phy <- gheatmap(phy, DF, low = 'white', high = 'blue', colnames = FALSE);
    if (!is.null(fname)) ggsave(plot = phy, filename = fname, ...);
    return(phy);
    }