cluster.nj <- function(dist.mat, n.cls) {
    nj.tree <- ape::bionj(as.dist(dist.mat));
    hc <- hclust(
        d = as.dist(cophenetic(nj.tree)),
        method = 'complete'
        );
    cls <- cutree(hc, k = n.cls);
    nj.tree$cluster <- cls[nj.tree$tip.label];
    nj.tree <- get.labels(tree = nj.tree);
    
    # measure average cluster 'tightness'
    nj.tree$asw <- summary(silhouette(cls, dist.mat), FUN = mean)$clus.avg.widths;

    # measure similarity between true clones and cluster
    nj.tree$ari <- mclust::adjustedRandIndex(nj.tree$true.label, nj.tree$cluster);
    return(nj.tree);
    }
