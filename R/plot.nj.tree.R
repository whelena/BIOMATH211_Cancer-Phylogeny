plot.nj.tree <- function(D.mat, fname = 'phylo_tree.png') {
    nj.tree <- ape::bionj(as.dist(D.mat));
    png(
        file = fname,
        type = 'cairo',
        width = 600,
        height = 400
        );
    plot.phylo(nj.tree);
    nodelabels(gsub('_.*', '', nj.tree$tip.label));
    dev.off();
    return(as.phylo(nj.tree));
    }
