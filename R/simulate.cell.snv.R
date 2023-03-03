simulate.cell.snv <- function(DT, n.cell) {
    n.snv <- length(unique(DT$snv.id));
    sim.cell <- setnames(data.table(
        matrix(apply(DT, 1, function(x) rbinom(n.cell, 1, as.numeric(x['CCF']))), n.cell, n.snv)
        ), DT$snv.id);
    }