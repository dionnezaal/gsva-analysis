## Bootstrapping gene sets
Gene sets are resampled with the same size as the initial gene sets. Genes from the original gene sets are used to fill the bootstrap gene sets, according to the same distribution as the genes in the original gene sets. When genes are used in a bootstrap gene sets they will be removed from the full list of genes, these will not be used again. Within the same gene sets no duplicate genes are allowed, in different gene sets the same genes can be used.


## P-value calculation
The p-value for the GSVA score indicates the certainty of the GSVA score. There is tested if it is possible to also get the same score when you would calculate a score for a random gene set of the same size. There is tested how often the score is better than when using random gene sets of the same size. 

The GSVA score calculated with original gene sets will be called the original GSVA score. The GSVA scores that are from the resampling of the gene sets will be called the bootstrap GSVA scores. 

The p-value is calculated by first separating the negative and positive scores from the bootstrapped scores. When the original score is positive (>=0) then there will be continued with only the positive scores, if the original score is negative (<0) then there will be continued with the negative scores. With either the positive or negative bootstrap scores there is assessed if the bootstrap scores is the same or higher than the original scores. The amount of times the bootstrap scores is higher or equal to the original score is divided by the amount of positive or negative score. This will result in a p-value score for the gene set. This p-value is already corrected for the amount of times the bootstrap is done, because there is divided by the amount of positive or negative scores there are.
