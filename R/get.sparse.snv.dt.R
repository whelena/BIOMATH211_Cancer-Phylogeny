get.sparse.snv.matrix <- function(DT) {
    ccf.array        <- as.data.frame(convert.df2array(DT));
    ccf.array$snv.id <- rownames(ccf.array);

    snv.dt <- melt(
        data          = setDT(ccf.array),
        id.vars       = c('snv.id'),
        measure.vars  = colnames(ccf.array)[which(colnames(ccf.array) != 'snv.id')],
        variable.name = 'ID',
        value.name    = 'CCF'
        );

    snv.dt <- merge(
        snv.dt,
        unique(DT[, !colnames(DT) %in% c('CCF', 'ID')]),
        by = c('snv.id'),
        all.x = TRUE,
        allow.cartesian = TRUE
        );
    return(unique(snv.dt));
    }
