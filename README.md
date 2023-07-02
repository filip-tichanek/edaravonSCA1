# edaravonSCA1
This page shows all data and full [R](https://www.r-project.org/) code for the study [Sucha et al.](https://www.mdpi.com/1422-0067/24/13/10689) (2023, *International Journal of Molecular Science*).

**Citation**: 
> Sucha M, Benediktova S, Tichanek F, Jedlicka J, Kapl S, Jelinkova D, Purkartova Z, Tuma J, Kuncova J, Cendelin J. Experimental Treatment with Edaravone in a Mouse Model of Spinocerebellar Ataxia 1. International Journal of Molecular Sciences. 2023; 24(13):10689. https://doi.org/10.3390/ijms241310689

This page shows HTML files with commented code and its output, and R code in quarto (.qmd) file. To explore the code and results, download the html file and open it in your browser, or click on link below (you will see the document on *RPubs* page). You can download it by clicking on the file and then clicking on 'download raw file' (icon for downloading on the right)

Bayesian regression was used with [brms](https://cran.r-project.org/web/packages/brms/index.html) R package

The following codes are shown:

- [***code01_data_functions***](https://rpubs.com/filip_tichanek/edaSCA1_code01) imports data, packages, and defines custom functions
- ***code02_behav*** shows an analysis of emotional-related behavioural indicators using  *Forced swim*, *sucrose preference* and *Open field* tests. The analysis of the 1st outcome (*Forced swim*) **provides a detailed guide** on how to interpret the outputs of Bayesian diagnostics, *posterior predictive check* as well as resulting posterior distributions. Thus, reading this will help to understand the next codes. Models were fitted either with Gaussian, gamma (with log-link, only-positive right-tailed outcomes) or beta (logit-link, continuous proportions/probabilities) likelihood distribution
- ***code03_grip_strength*** shows analysis of *grip strength* with gamma model with log-link and with the factor of the individual mouse (*id*) as a random intercept
- ***code04_CatWalk*** shows an analysis of Gait Abnormalities. Besides Bayesian regression (Gaussian multivariate regression), this analysis was simultaneously done in a classical (frequentist) framework using PERMANOVA. 
- ***code05_rotarod*** shows an analysis of motor coordination on the *Rotarod* test using hierarchical (mixed-effects) Bayesian Gamma (log-linked) regression
- ***code06_T_maze*** shows analysis of cognition (learning and cognitive flexibility) in *T-maze*. The summarized results are analysed as univariate outcomes via beta-binomial regression
- ***code07_elisa*** analysis of hippocampal and cerebellar levels of interleukin 6 (IL6) and cerebellar levels of brain-derived neurotrophic factor (BDNF). Analysed via Gaussian regression
- ***code08_mitochondria_respiration*** shows an analysis of mitochondrial respiration capacity in both the cerebellar (*cb*) and hippocampal (*hp*) tissue. Respiration was standardized by either tissue wet weight of (*mg*) or by citrate synthase activity (*cs*), across different phases experiment. Data per tissue and standardization were analysed separately, using multivariate Gamma models with log-link, where respiration in different states of the substrate-uncoupler-inhibitor protocol represented multiple outcomes.
- ***code09_citrate_syn*** analyze citrate synthase activity in both cerebellar and hippocampal tissue, using gamma models with log-link 
- ***code10_calbindin*** shows an analysis of calbindin immunofluorescence intensities in the cerebellar molecular layer. The analysis was done via a hierarchical Gaussian model (with a random intercept of *slice* nested in the random intercept of *id*/*subject*), and with the cerebellar region (hemisphere vs. vermis) and immunofluorescence intensity in the neighbouring granular layer as covariates. 
- ***code11_cb_volume*** analyze cerebellar molecular layer volume. As data exploration suggested outlying values (of both directions), robust regression with student t-distribution and with random intercept (*id*) was applied.
- ***code12_hp_volume*** analyzes hippocampal volumes via robust regression with student t-distributed likelihood and with random-intercept of *id*
- ***code13_psa_ncam*** shows analysis of PSA-NCAM immunofluorescence intensity, presumably reflecting processes related to neuroplasticity. We used Bayesian hierarchical generalized additive model with Gamma distribution and random intercept (random effect: mouse *id*). As there was an apparent sex-related difference, further confirmed with leave-one-out cross-validation, the sex factor was included in the final model. As the PSA-NCAM IF varied across slices from the frontal to caudal regions, the order of slices was included as a factor with a non-linear effect, fitted with thin-plate splines limited to 3 knots.

From code02 and further, the codes have the following structure:

- general section **Import of packages, data and functions**, loading file saved by *code01_data_functions*, uploading packages used and providing some additional data wrangling (modification data to the form suitable for downstream analysis). 

  The next sections define the specific outcome to be modelled (e.g., distance walked in an open field over 10 minutes) and include:
  - **data exploration** showing data distribution per group, and searching for possible covariates that need to be adjusted (e.g. sex)
  - **modelling** with the sections
      - priors specification
      - model(s) fitting
      - model diagnostics (searching for signs of problems with convergence or insufficient sample size)
      - [posterior predictive check](https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html) to explore if simulations from the posterior distributions reconstruct our data well (whether the model fits data well and the distributional assumptions are met).
      - extraction of posterior samples (needed for visualisation)
  - **visualisation** with code showing
      - data and group-specific predictions (and their 95% credible interval)
      - posterior distribution of effect size (between-groups differences)
    

### Priors specification

Prior probability distributions for all fixed-effect parameters of Bayesian models were explicitly specified to have normal distribution on the natural scale of a given model (log for gamma models, logit for [beta-]binomial and beta models).

Priors for the effect of edaravone treatment were set to have a mean value equal to zero (*mu = 0*), thus attributing the highest prior probability to the absence of an effect. In log- and logit-linked models, the prior *sigma* was set to 1.2. In Gaussian and robust models, prior sigma was set to 1.2 of standard deviations of the outcome in SCA1-specific data subset (sigma = 1.2 x SD[outcome *i*]), addressing the expectancy that the edaravone may improve but highly likely not resolve the SCA1 phenotype. 

Priors for the effect of genotype were derived from our previous results of behavioral and neuropathological abnormalities in SCA1 mice [Tichanek et al., 2020, Scientific Reports](https://www.nature.com/articles/s41598-020-62308-0) when there was “statistically significant” genotype-related difference in the previous data. *Mu* was set to be equal to the estimate based on the previous data (β<sub>2020</sub>) ± its standard error (SE[<sub>2020</sub>]; *mu* = β<sub>2020</sub>  ±  SE[β<sub>2020</sub>] in the direction toward the zero effect. Sigma was set to be 1.5-fold of the distance of the mu to either β<sub>2020</sub> or zero (the higher from these two distances). Only results from the mice younger than 16 weeks of age were used for prior distribution specification. 

For analyses which have not been performed in the previous manuscript, or when methodologies were not the same, we used  *mu = 0* and *sigma* of 2 x SD [outcome *i*] (SD of the outcome *i* in either control or SCA1 mice, depending on which SD was higher).
Priors for intercept (control SCA1 mice were set as a reference group in all models) were based to be between the previous SCA1 mean from the above-cited data and the overall mean from current data. Sigma was set to 5 x SD [outcome *i*]


### General concept of the study 

The original study aimed to explore the potential of edaravone (mitochondrial-targeted drug) to slow down the progression of the neurodegenerative disease [Spinocerebellar Ataxia 1](https://en.wikipedia.org/wiki/Spinocerebellar_ataxia_type_1) (SCA1) in mice genetically engineered (*knock-in*) to show hallmarks of the disease (***SCA1 mice***)

Thus, the main objective of all analyses is the comparison of edaravone-treated vs. control (on saline) mice with SCA1 disease in numerous indicators of the severity of the disease (although data on the effects of the presence of the disease and the impact of edaravone on healthy mice were also analysed for most outcomes). In sum, we had the following 3 (4) groups of mice

- healthy mice treated with saline (***WT_0***)
- healthy mice treated with edaravone (***WT_E***, *! only for behavioural and mitochondrial-function-related outcomes*)
- SCA1 mice treated with saline (***SCA1_0***)
- SCA1 mice treated with edaravone (***SCA1_E***)
