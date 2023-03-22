library(ape);
library(data.table);

name  <- 'sim_K3_S3_T1000_M30_G0_run3';
repo.dir <- '/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny'
sc.dir <- file.path(repo.dir, 'data', 'input', 'simulated_single_cell');
out.dir <- file.path(repo.dir, 'data', 'output', 'name');
devtools::load_all(repo.dir);

sim.sc <- fread(
    file = file.path(sc.dir, paste(name, N.total, 'sc.tsv', sep = '-')),
    sep = '\t',
    header = TRUE
    );


tree <- plot.nj.tree(
    D.mat = dist(sim.sample.dt[, 1:100], method = 'manhattan'),
    fname = file.path(out.dir, 'sc_tree.png')
    );