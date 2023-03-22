simulate.cell.snv <- function(DT, n.cell, sample.name = NULL) {
    n.snv <- length(unique(DT$snv.id));
    sim.cell <- setnames(data.table(
        matrix(apply(DT, 1, function(x) sort(rbinom(n.cell, 1, as.numeric(x['CCF'])))), n.cell, n.snv)
        ), unique(DT$snv.id));
    if (!is.null(sample.name)) {
        sim.cell[, sample.id := sample.name];
        }
    return(sim.cell)
    }
# improvement is to use a zero inflated binomial distribution to set a non-zero probability that a cell has a  mutation even when the given VAF is zero

simulate.cell.per.sample <- function(DT, n.cell) {
    n.sample <- length(unique(DT$ID));
    sim.list <- list();
    for (s in unique(DT$ID)) {
        sim.list[[s]] <- simulate.cell.snv(
            DT = copy(DT)[ID == s, ],
            n.cell = ceiling(n.cell / n.sample),
            sample.name = s
            );
        }
    sim.dt <- do.call('rbind', sim.list);
    }

simulate.single.cell <- function(
    n.cell = 10,
    n.clone = 3,
    n.snv = 10,
    noise = 0.05
    ) {
    
    clone.cp <- simulate.nested(
        n = n.clone,
        max.val = 1,
        min.val = 0,
        name.prefix = 'clone'
        );
    clone.snvs <- simulate.nested(
        n = n.clone,
        max.val = n.snv,
        min.val = 1,
        name.prefix = 'clone',
        diff = FALSE,
        integer = TRUE
        );
    
    cells <- list();
    for (i in seq_along(clone.snvs)) {
        for (j in seq_len(n.cell * clone.cp[i])) {
            cells[[paste(names(clone.snvs)[i], j , sep = '_')]] <- c(
                rbinom(clone.snvs[i], 1, (1 - noise)),
                rbinom((n.snv - clone.snvs[i]), 1, noise / 2)
                );
            }
        }
    cell.dt <- data.table(do.call('rbind', cells));
    setnames(cell.dt, paste0('snv.', 1:length(cell.dt)));
    setkey(cell.dt[, cell.id := names(cells)], cell.id);

    clone <- data.table(
        parent.id = 0:(n.clone - 1),
        clone.id = 1:n.clone,
        CP = rev(cumsum(rev(clone.cp))),
        distance = c(clone.snvs[1], diff(clone.snvs))
        );
    return(list('cells' = cell.dt, 'clone' = clone));
    }

simulate.nested <- function(
    n = 3,
    max.val = 1,
    min.val = 0,
    name.prefix = 'clone',
    diff = TRUE,
    integer = FALSE
    ) {

    sim.list <- c(max.val);
    for (i in 2:n) {
        max.val <- runif(1, min.val, max.val);
        if (integer) {
            max.val <- floor(max.val);
            }
        sim.list <- c(max.val, sim.list);
        }

    if (diff) {
        sim.list <- c(max.val, diff(sim.list))
        }
    names(sim.list) <- paste(name.prefix, 1:n, sep = '.');
    return(sim.list);
    }
