# overlapper

A shiny app to explore overlaps between lists. Includes:

* Venn diagrams
* UpSet plots
* Dotplots of pairwise overlap (with enrichment, size and significance of the overlaps)
* An interactive heatmap of pairwise overlaps reporting in addition the actual items of the intersections
* Permutation-based significance tests of the overlap between more than 2 sets.

Install with
```
BiocManager::install("plger/overlapper")
```

and then launch with
```
library(overlapper)
overlapper()
```
