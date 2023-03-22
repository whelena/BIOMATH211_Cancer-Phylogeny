library(ape);
library(ggtree);
library(argparse);
library(RColorBrewer);
repo.dir <- '/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny';
devtools::load_all(repo.dir);
devtools::load_all('/hot/user/hwinata/MSK_cfDNA/public-R-CancerEvolutionVisualization/');

parser <- ArgumentParser();
parser$add_argument(
    '-c', '--nclone',
    help = 'name of simulated data'
    );
parser$add_argument(
    '-s', '--nsnv',
    help = 'name of simulated data'
    );
parser$add_argument(
    '-n', '--ncell',
    help = 'name of simulated data'
    );
parser$add_argument(
    '-b', '--noise',
    default = 0.95
    );
ARGS <- parser$parse_args()

# ARGS <- list(
#     nclone = 3,
#     nsnv = 100,
#     ncell = 100,
#     noise = 0.3
#     );
for (var in names(ARGS)) {
    assign(var, as.numeric(ARGS[[var]]));
    }
name <- paste0('sim_K', nclone, '_M', nsnv, '_C', ncell, '_N', noise);
out.dir <- file.path(repo.dir, 'data',  name);
if (!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE);
###################################################################################################
# Simulate single cell calls
###################################################################################################
sim <- simulate.single.cell(
    n.cell = ncell,
    n.clone = nclone,
    n.snv = nsnv,
    noise = noise
    );

cell.dt <- sim$cells;
ancs.dt <- sim$clone;

fwrite(
    x = cell.dt,
    file = file.path(out.dir, 'simulated_sc.tsv'),
    sep = '\t',
    quote = FALSE,
    col.names = TRUE
    );

###################################################################################################
# NJ - clustering and visualization
###################################################################################################
dist.mat <- proxy::dist(
    x = cell.dt[, paste0('snv.', 1:nsnv), with = FALSE],
    method = 'Jaccard'
    );
dimnames(dist.mat) <- cell.dt$cell.id;
if (any(dist.mat > 100)) {
    dist.mat <- (scale(dist.mat) + 1) * 25;    
    }
nj.tree <- cluster.nj(dist.mat, n.cls = nclone);
plot.figure(
    tree = nj.tree,
    DT = cell.dt,
    fname = file.path(out.dir, 'tree.png')
    );
###################################################################################################
# Save stats
###################################################################################################
cell.dt[, `:=`(
    true.clone.id = gsub('_.*', '', cell.id),
    clone.id = nj.tree$cluster
    )];
ancs.dt[, `:=`(
    sim.cp = rev(cumsum(rev(cell.dt[, .N / ncell, by = true.clone.id]$V1))),
    cls.cp = rev(cumsum(rev(cell.dt[, .N / ncell, by = clone.id]$V1))),
    asw = nj.tree$asw
    )];
fwrite(
    x = ancs.dt,
    file = file.path(out.dir, 'simulated_phylogeny.tsv'),
    sep = '\t',
    quote = FALSE,
    col.names = TRUE
    );
print(paste('ARI:', nj.tree$ari));

########################################################################################################
# Plot Tree
########################################################################################################
tree.dt <- data.table(
    parent  = ancs.dt$parent.id,
    label   = ancs.dt$clone.id,
    length1 = ancs.dt$distance
    );
plt <- SRCGrob(
    tree.dt,
    main = paste(name, 'Tree'),
    main.cex = 1,
    node.col = unique(nj.tree$true.col),
    horizontal.padding = 0,
    scale1 = 2,
    yaxis1.label = 'Distance (# of SNV)',
    add.normal = TRUE
    );

png(
    file = file.path(out.dir, 'CEV-tree.png'),
    width = 5,
    height = 7,
    units = 'in',
    res = 200
    );
grid.draw(plt);
dev.off();
