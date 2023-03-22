find.ancestor <- function(child.id, phy.dt) {
    # Find the parent(s) of the current child
    parents <- phy.dt[phy.dt$child == child.id, "parent"];

    # If there are no parents, return an empty list
    if (0 == parents | is.na(parents)) {
        return(list());
        }

    # Otherwise, recursively find the ancestors of each parent and combine the results
    ancestor.list <- list()
    for (parent.id in parents) {
        parent.ancestor <- find.ancestor(parent.id, phy.dt)
        ancestor.list <- c(ancestor.list, parent.id, parent.ancestor)
        }
    return(unique(ancestor.list));
    }
