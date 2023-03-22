get.labels <- function(tree) {
    tree$true.label <- gsub('_.*', '', tree$tip.label);
    tree$true.col <- get.colours(value.list = tree$true.label);
    return(tree);
    };


get.colours <- function(
    value.list,
    return.names = FALSE
    ) {

    n        <- length(unique(value.list));
    col.list <- c(brewer.pal(8, 'Dark2'), brewer.pal(5, 'Set1'))[1:n];
    if (is.null(levels(value.list))) {
        value.list <- factor(value.list, levels = unique(value.list))
        }
    names(col.list) <- levels(value.list);
    if (return.names) {
        return(col.list);
    } else {
        return(col.list[value.list]);
        }
    }
