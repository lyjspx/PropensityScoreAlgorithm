# Propensity Score Algorithm

### What is Propensity Score Algorithm?

The propensity score algorithm is a statistical method used in observational studies to control for potential confounding variables. In observational studies, researchers cannot randomly assign subjects to different groups, so there may be differences in the baseline characteristics of the groups being compared. The propensity score algorithm attempts to balance these differences by creating a single score that represents the probability of being assigned to a particular group based on a set of pre-treatment variables.

The basic steps of the propensity score algorithm are as follows:

1. Determine the pre-treatment variables that may be associated with the outcome of interest.
2. Use a logistic regression model to estimate the probability of being assigned to the treatment group (i.e., the propensity score) based on the pre-treatment variables.
3. Check the balance of the pre-treatment variables between the treatment and control groups using standardized differences or other statistical tests.
4. Use the propensity score to adjust for imbalances in the pre-treatment variables by matching or weighting the treatment and control groups based on their propensity score.
5. Assess the balance of the pre-treatment variables after adjustment to ensure that they are now comparable between the treatment and control groups.

The propensity score algorithm can be a useful tool to reduce the bias in observational studies and to estimate the effect of a treatment or intervention. However, it is important to carefully select the pre-treatment variables and to check the balance of these variables before and after adjustment to ensure that the results are valid.

### What about this repository?

The original script, sourced from https://github.com/youqiongye/SexImm, has limited transferability to other analyses due to some out-of-dated packages and its non-generic characteristics. To address this, the repository seeks to update the packages and provide generic interfaces for performing Propensity Score Algorithm.

Citations:

Ye Y, Hu Q, Chen H, Liang K, Yuan Y, Xiang Y, Ruan H, Zhang Z, Song A, Zhang H, Liu L, Diao L, Lou Y, Zhou B, Wang L, Zhou S, Gao J, Jonasch E, Lin SH, Xia Y, Lin C, Yang L, Mills GB, Liang H, Han L. Characterization of Hypoxia-associated Molecular Features to Aid Hypoxia-Targeted Therapy. ***Nat Metab.*** 2019 Apr;1(4):431-444. doi: 10.1038/s42255-019-0045-8. Epub 2019 Mar 18. PMID: 31984309; PMCID: PMC6980239.

Li, L. and Greene, T., 2013. A weighting analogue to pair matching in propensity score analysis. ***The international journal of biostatistics***, *9*(2), pp.215-234.