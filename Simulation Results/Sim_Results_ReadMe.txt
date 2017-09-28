This ReadMe file provides a brief summary of the results found in the given text files.

eb10_constant, eb20_constant, equal_constant, pooled_constant:  these files are the results from gradient descent style algorithm to identify the ideal thresholds for the approaches incorporating supplemental information when we are attempting to control the type-1 error rate under the constant mortality scenario. The final row indicates the thresholds to be used for each segment.

simresults_mmddyy:  contains the xtable formatted output included in the manuscript and supplementary materials, however further formatting was done for publication by combining different parts of the tables together (i.e., these files show constant and varying mortality results in the same table for each scenario, but we primarily present them as all constant mortality results in one table (e.g., Table 2 in the manuscript) and all varying mortality results in one table (e.g, Table 3 in the manuscript)

simresults_futility_mmddyy:  similar to simresults_mmddyy, but for scenarios which incorporate a posterior probability futility boundary

simresults_constant, simresults_varying:  the raw output from simulations, primarily of use to create the barplot figures in the "Figures_Code.R" file