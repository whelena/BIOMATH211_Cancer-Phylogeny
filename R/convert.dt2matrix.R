convert.dt2matrix<- function(
    DT,
    value = 'CCF',
    x.axis = 'snv.id',
    y.axis = 'ID'
    ) {

    if (is.null(DT[[x.axis]]) | is.null(DT[[value]]) | is.null(DT[[y.axis]])) {
        stop(paste('Dataframe does not contain one of the columns:', value, x.axis, y.axis));
        }

    DT <- setnames(
        x = DT[, c(x.axis, y.axis, value), with = FALSE],
        old = c(x.axis, y.axis, value),
        new = c('x', 'y', 'val')
        );
    DT  <- unique(DT);
    arr <- dcast(DT, x ~ y, value.var = 'val', fill = 0);
    # set x.axis as rownames
    rows            <- arr[[1]];
    cols            <- names(arr)[-1];
    arr             <- as.matrix(arr[, -1]);
    rownames(arr)   <- rows;
    colnames(arr)   <- cols;

    if (!is.null(levels(DT$y)) & ncol(arr) > 1) {
        arr <- arr[, levels(DT[, y])];
        }

    return(arr);
    }
