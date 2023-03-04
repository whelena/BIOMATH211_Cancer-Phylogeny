library(ape);
library(phangorn);
library(phytools);
library(data.table);

repo.dir <- '/Users/hwinata/Library/CloudStorage/OneDrive-UCLAITServices/2022-2023/Winter/BIOMATH 211/BIOMATH211_Cancer-Phylogeny'
out.dir <- file.path(repo.dir, 'data', 'output');
devtools::load_all(repo.dir);

# read data
snv.dt <- fread(
    file = 'data/input/K1_S3_M30.tsv',
    sep = '\t',
    header = TRUE
    );
snv.dt[, median.ccf := CCF, by = 'clone.id'];
# convert to taxa by seq matrix where taxa are clones and sequences are SNVs with
snv.mat <- convert.dt2matrix(
    DT = snv.dt,
    value = 'median.ccf',
    x.axis = 'clone.id',
    y.axis = 'snv.id'
    );

# sim.cell <- unique(sim.cell[, cell.count := .N, by = names(sim.cell)]);

sim.cell <- simulate.cell.snv(snv.dt, 1000);
tree1 <- plot.nj.tree(dist(unique(sim.cell)));

D.mat <- dist(sim.cell, method = 'euclidean');

cls <- cluster.cells(
    D.mat,
    linkage = 'complete',
    n.cls = 4,
    fname = file.path(out.dir, 'test-cls.png')
    )

# to get CP
table(cls) / 1000; 

cell.dt <- setDT(sim.cell)[, 'cls' := cls][, lapply(.SD, mean), by = 'cls', .SDcols = names(sim.cell)[1:30]];
tree <- plot.nj.tree(dist(cell.dt));
mrca(tree, full = TRUE);

parsimony(tree, as.phyDat(cell.dt[1:30]));
tree2 <- root(tree, out = 1)
plot(tree2, show.tip = TRUE, edge.width=2)

myBoots <- boot.phylo(tree, as.matrix(cell.dt[, 2:31]), function(e) root(bionj(dist(e)), 1), trees = TRUE);


########################################################################################################
# Migration pattern
snv.dt <- fread(
    file = 'data/input/patient7.tsv',
    sep = '\t',
    header = TRUE
    );

cls <- snv.dt$clone.id[!(grepl('^S..', snv.dt$clone.id) | is.na(snv.dt$clone.id))];
snv.dt <- snv.dt[clone.id %in% cls, ];

snv.mat <- convert.dt2matrix(
    DT = snv.dt,
    value = 'CCF',
    x.axis = 'ID',
    y.axis = 'snv.id'
    );

bin.snv.mat <- snv.mat;
bin.snv.mat[bin.snv.mat > 0] <- 1;
tree <- plot.nj.tree(dist(bin.snv.mat, method = 'manhattan'));
root(tree, out = c(1:3));

myBoots <- boot.phylo(tree, as.matrix(cell.dt[, 2:31]), function(e) root(bionj(dist(e)), 1), trees = TRUE);

