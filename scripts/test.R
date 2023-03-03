library(ape);
library(phangorn);
library(phytools);
library(data.table);

devtools::load_all('/Users/hwinata/Library/CloudStorage/OneDrive-UCLAITServices/2022-2023/Winter/BIOMATH 211/BIOMATH211_Cancer-Phylogeny');

# read data
snv.dt <- fread(
    file = 'data/K1_S3_M30.tsv',
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

CP <- apply(snv.mat, 1, function(x) mean(x[x > 0]));

snv.mat <- snv.mat[order(CP, decreasing = T), ];

combn(rownames(snv.mat), 2, function(x) length(setdiff(snv.mat[x[[1]], ], snv.mat[x[[2]], ])))

nj.tree <- ape::nj(as.dist(n2.mat));

plot.phylo(nj.tree);
nodelabels();
tiplabels();

N.cell <- 10;
N.snv <- length(unique(snv.dt$snv.id));
sim.cell <- setnames(data.table(
    matrix(apply(snv.dt, 1, function(x) rbinom(N.cell, 1, as.numeric(x['CCF']))), N.cell, N.snv)
    ), snv.dt$snv.id);

sim.cell <- unique(sim.cell[, cell.count := .N, by = names(sim.cell)]);

nj.tree <- ape::nj(dist(sim.cell[, 1:N.snv]));
plot.phylo(nj.tree);
nodelabels();

simulate.cell.snv <- function(DT, n.cell) {
    n.snv <- length(unique(snv.dt$snv.id));

}