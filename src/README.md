## bug fixes

* The default output directory in base ReVis -> currently the default if it is not specified is supposed to be the current working directory, but that somehow adds "None" as a prefix in front of the filename? haven't figured out yet where that shows up.

## improvement suggestions

### Base ReVis

* Add an option to give custom names instead of contig names to the contigs in the plot

### Gene surroundings ReVis

* !! implement a runtimeError if the .ori.out (or .out) file contains contigs not present in the assembly! currently only counts 0 repeats if that is the case.
  * In my *C. maculatus* not superscaffolded repeatmasked assembly (ENA download) in .ori.out the contig names are somehow parsed wrong and all are called `ENA` instead of the actual `utg000[...]` name in the assembly and annotation

**Assess enrichment through CI overlap:** [This blog post](https://statisticsbyjim.com/hypothesis-testing/confidence-intervals-compare-means/) explains it quite nicely, but I should make the confidence interval 83% instead of 95%, see [Goldstein and Healy (1995)](https://rss.onlinelibrary.wiley.com/doi/abs/10.2307/2983411).