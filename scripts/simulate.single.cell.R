library(ape);
library(ggtree);
library(argparse);
repo.dir <- '/hot/user/hwinata/BIOMATH211_Cancer-Phylogeny';
devtools::load_all(repo.dir);


parser <- ArgumentParser();
parser$add_argument(
    '-n', '--name',
    help = 'name of simulated data'
    );
name <- parser$parse_args()$name;

name  <- 'sim_K3_S10_T1000_M30_G0_run3';
sim.dir <- file.path(repo.dir, 'data', 'input', 'simulation');
sc.dir <- file.path(repo.dir, 'data', 'input', 'simulated_single_cell');
out.dir <- file.path(repo.dir, 'data', 'output', name);
if (!file.exists(sc.dir)) dir.create(sc.dir, recursive = TRUE);
if (!file.exists(out.dir)) dir.create(out.dir, recursive = TRUE);

# read data: 10 clones, 3 samples, 100 snvs
snv.dt <- fread(
    file = file.path(sim.dir, paste0(name, '-src_input.tsv')),
    sep = '\t',
    header = TRUE
    );
snvs <- unique(snv.dt$snv.id);
tree.dt <- fread(
    file = file.path(sim.dir, paste0(name, '-structure.tsv')),
    sep = '\t',
    header = TRUE
    );

leafs <- setdiff(tree.dt$child, tree.dt$parent);
ancestor.dt <- data.table(child = tree.dt$child);
ancestor.dt$ancestors <- lapply(ancestor.dt$child, find.ancestor, phy.dt = tree.dt);
###################################################################################################
# Simulate single cell calls
###################################################################################################
N.total <- 100;
n.snv <- length(unique(snv.dt$snv.id));
sim.list <- list();
for (i in leafs) {
    cls <- unname(unlist(ancestor.dt[child == i,]));
    tmp <- copy(snv.dt)[!clone.id %in% cls, CCF := 0];
    sim.list[[i]] <- simulate.cell.per.sample(
        DT = tmp,
        n.cell = ceiling(N.total / length(leafs))
        );
    }

sim.dt <- do.call('rbind', sim.list);
sim.dt <- sim.dt[rowSums(sim.dt[, ..snvs]) > 0, ];
setkey(sim.dt[, cell.id := paste0('cell.', 1:nrow(sim.dt))], cell.id);

fwrite(
    x = sim.dt,
    file = file.path(sc.dir, paste(name, N.total, 'sc.tsv', sep = '-')),
    sep = '\t',
    quote = FALSE,
    col.names = TRUE
    ); 
###################################################################################################
# get clone mutation profile
###################################################################################################
clone.dt <- dcast(
    data = snv.dt,
    formula = snv.id ~ clone.id,
    value.var = 'CCF',
    fill = 0,
    fun.aggregate = max
    );
ancestor.dt[, c(snvs) := 0];
for (i in ancestor.dt$child) {
    cls <- unique(as.character(unname(unlist(ancestor.dt[child == i, c(child, ancestors)]))));
    ancestor.dt[child == i, c(snvs) := as.list(rowSums(clone.dt[, cls, with = FALSE]) / length(cls))];
    }
ancestor.dt[, c(snvs) := ceiling(.SD), .SDcols = c(snvs)];

fwrite(
    x = ancestor.dt[, ancestors := NULL],
    file = file.path(sc.dir, paste0(name, '-clone_mut_profile.tsv')),
    sep = '\t',
    quote = FALSE,
    col.names = TRUE
    );

###################################################################################################
# get similarity measure between each cell to each clone
###################################################################################################
score.dt <- sim.dt[, c('sample.id', 'cell.id')];
for (i in score.dt$cell.id) {
    cell.mp <- as.matrix(sim.dt[cell.id == i, ..snvs])[1, ];
    clone.mp <- as.matrix(ancestor.dt[, ..snvs]);
    score <- (n.snv - rowSums(abs(cell.mp - clone.mp))) / n.snv;
    score.dt[cell.id == i, paste('clone', ancestor.dt$child, sep = '.') := as.list(score)];
    }
score.dt[, `:=`(
    min.dist.clone = names(.SD)[max.col(.SD)],
    score = do.call(pmax, .SD)
    ), .SDcols = paste('clone', ancestor.dt$child, sep = '.')
    ];
fwrite(
    x = score.dt,
    file = file.path(out.dir, 'score.tsv'),
    sep = '\t',
    quote = FALSE,
    col.names = TRUE
    );
###################################################################################################
# PLOTTING
###################################################################################################
score.dt[, `:=`(
    N = .N,
    mean.score = mean(score)
    ), by = .(sample.id, min.dist.clone)
    ];
bar.dt <- unique(score.dt[, c('sample.id', 'min.dist.clone', 'N', 'mean.score')])

phy.dt <- rbind(
    sim.dt[, ID := cell.id][, c('ID', snvs), with = FALSE],
    ancestor.dt[, ID := paste0('clone.', child)][, c('ID', snvs), with = FALSE]
    );

cell.dist <- dist(phy.dt[, ..snvs], method = 'euclidean');
dimnames(cell.dist) <- phy.dt$ID;

nj <- plot.nj.tree(
    D.mat = cell.dist,
    fname = file.path(out.dir, 'nj-tree.png')
    );

nj$tip.color <- unlist(lapply(nj$tip.label, function(x) if (grepl('cell', x)) 'steelblue' else 'salmon'));
phy <- ggtree(as.phylo(nj)) + geom_tippoint(color = nj$tip.color, size = 3);
 geom_facet(panel = 'clone', data = score.dt[, c('cell.id', 'min.dist.clone')], geom = geom_point, 
        mapping = aes(x = 1, color = min.dist.clone), shape = '|', size = 5)
# phy.long <- melt(
#     data = phy.dt,
#     id.vars = 'ID',
#     measure = snvs,
#     variable.name = 'SNV',
#     value.name = 'presence'
#     );
# phy.long[, SNV := as.numeric(gsub('s', '', SNV))];
# phy + geom_facet(panel = "SNP", data = phy.long, geom_point, mapping = aes(x = SNV, color = presence), shape = '|')


###
n.cells <- 100;
clone.prop <- c(0.1, 0.2, 0.4, 0.3);
ccfs <- cumsum(clone.prop);
names(clone.prop) <- c('a', 'b', 'c', 'd');
snv.max <- 100;
snv.list <- list();
for (clone in names(clone.prop)) {
    snv.list[clone] <- floor(runif(1, 1, snv.max));
    snv.max <- snv.list[[clone]];
    }
clone.snvs <- unlist(snv.list);
names(clone.snvs) <- names(clone.prop);
cells <- list();
for (j in seq_along(clone.snvs)) {
    for (i in seq_len(n.cells * clone.prop[j])) {
        cell <- rbinom(clone.snvs[j], 1, 0.95);
        if (j > 1) cell <- c(rep(0, (clone.snvs[1] - clone.snvs[j])), cell);
        cells[[paste(names(clone.prop)[j], i, sep = '.')]] <- cell;
        }
    }
cells.dt <- as.data.frame(t(as.data.frame(cells)));
ccfs.simulated <- sapply(cells.dt, mean);
cell.dist <- dist(cells.dt, method = 'euclidean');
dimnames(cell.dist) <- rownames(cells.dt)
nj <- plot.nj.tree(
    D.mat = cell.dist,
    fname = file.path(out.dir, 'nj-tree.png')
    );
clone.label <- gsub('\\..*', '', rownames(cells.dt));
col.scheme <- list('a' = 'steelblue', 'b' = 'forestgreen', 'c' = 'maroon', 'd' = 'purple');
nj$tip.color <- unlist(col.scheme[clone.label]);
phy <- ggtree(as.phylo(nj)) + geom_tippoint(color = nj$tip.color, size = 3);

gheatmap(phy, cells.dt, low = 'white', high = 'blue', colnames = FALSE)
###
cutree(chronos(nj), k = 3)

hm <- as.matrix(phy.dt[, ..snvs]);
rownames(hm) <- phy.dt$ID;

hm2 <- data.frame(Clone = c(score.dt$min.dist.clone, ancestor.dt$ID));
rownames(hm2) <- c(score.dt$cell.id, ancestor.dt$ID)
phy <- gheatmap(phy, hm, low = 'white', high = 'blue', colnames = FALSE);
# phy <- gheatmap(phy, hm2, width = 0.1) +  scale_fill_manual(breaks = ancestor.dt$ID, 
#         values = c("steelblue", "firebrick", "darkgreen"), name = "Clone")
score.long <- melt(
    data = score.dt,
    id.vars = 'cell.id',
    measure = ancestor.dt$ID,
    variable.name = 'clone',
    value.name = 'score'
    );
phy <- phy + geom_facet(panel = 'Score', data = score.long, geom = geom_col,
    mapping = aes(x = score, color = as.factor(clone)), orientation = 'y');

phy <- phy + geom_facet(panel = 'clone', data = score.dt[, c('cell.id', 'min.dist.clone')], geom = geom_point, 
        mapping=aes(x = 1, color = min.dist.clone), shape = '|', size = 5)

ggsave(plot = phy, filename= file.path(out.dir, 'nj-ggtree.png'));
