plot.nj.tree <- function(D.mat) {
    nj.tree <- ape::bionj(as.dist(D.mat));
    plot.phylo(nj.tree);
    return(as.phylo(nj.tree));
    }
