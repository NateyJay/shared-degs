A function to classify genes from a set of DE-seq analyses as "shared" or not. This uses an empirical-FDR to generate a threshold of experiments. Genes that are DE in a given direction for this threshold number are considered a "shared DEG".

The function takes a simple input: a table of LFC values. 
* **Columns** represent experiments, which should follow a naming convention including what study the data is derived from (to prevent over-influence of libraries from the same study). This format is "Study.experiment_ID".
* **Rows** represent genes by name.
* The data frame can only contain numeric values (and NAs).

This function returns a list, which contains general information on the analysis, aggregated counts of DE genes in both directions, and a summary table which includes all of the shared deg values, including a call for which direction they are DE.

