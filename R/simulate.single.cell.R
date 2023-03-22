simulate.single.cell <- function(
    n.cell = 10,
    n.clone = 3,
    n.snv = 10,
    noise = 0.95
    ) {
    
    clone <- simulate.clone(n.clone, n.snv, noise);
    cell.dt <- simulate.cell(clone.cp = clone$cp, clone.snv = clone$snv, n.cell, n.snv);

    clone.dt <- data.table(
        parent.id = 0:(n.clone - 1),
        clone.id = 1:n.clone,
        CP = rev(cumsum(rev(clone$cp))),
        distance = c(clone$snv[1], diff(clone$snv))
        );
    return(list('cells' = cell.dt, 'clone' = clone.dt));
    }

simulate.clone <- function(
    n.clone = 3,
    n.snv = 10,
    noise = 0.95
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
    return(list('cp' = clone.cp, 'snv' = clone.snvs));
    }

simulate.cell <- function(
    clone.cp,
    clone.snv,
    n.cell,
    n.snv
    ) {
    cells <- list();
    for (i in seq_along(clone.snv)) {
        for (j in seq_len(n.cell * clone.cp[i])) {
            cells[[paste(names(clone.snv)[i], j , sep = '_')]] <- c(
                rbinom(clone.snv[i], 1, (1 - noise)),
                rbinom((n.snv - clone.snv[i]), 1, (noise / 2))
                );
            }
        }
    cell.dt <- data.table(do.call('rbind', cells));
    setnames(cell.dt, paste0('snv.', 1:length(cell.dt)));
    setkey(cell.dt[, cell.id := names(cells)], cell.id);
    return(cell.dt)
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
        max.val <- runif(n = 1, min = min.val, max = max.val);
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
