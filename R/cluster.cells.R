cluster.cells <- function(
    D.mat,
    linkage = 'complete',
    n.cls = 3,
    fname = 'cls.png'
    ) {

    hc  <- hclust(as.dist(D.mat), method = linkage);
    cls <- cutree(tree = hc, k = n.cls);
    png(
        file = fname,
        type = 'cairo',
        width = 600,
        height = 400
        );
    plot(
        x = hc,
        labels = row.names(hc),
        cex = 0.5,
        hang = 0
        );
    rect.hclust(tree = hc, k = n.cls, cluster = cls);
    dev.off();
    return(cls);
    }
