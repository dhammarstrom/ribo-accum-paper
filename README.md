# Ribosome accumulation during early phase resistance training

Daniel Hammarström<sup>1,2,£</sup>&ast;,Sjur J. Øfsteng<sup>2,£</sup>, Nicolai B. Jacobsen<sup>1</sup>, Krister B. Flobergseter<sup>1</sup>,  Bent R. Rønnestad<sup>1</sup>, Stian Ellefsen<sup>1,3</sup>


<sup>1</sup> Section for Health and Exercise Physiology, Department of Public Health and Sport Sciences, Inland Norway University of Applied Sciences, Lillehammer, Norway. 
<sup>2</sup> Swedish School of Sport and Health Sciences, Stockholm, Sweden.
<sup>3</sup> Innlandet Hospital Trust, Lillehammer, Norway. </div>
<sup>£</sup> These authors contributed equally to this work

----------------------------------------

This repository contains data, code and additional files needed to reproduce the analyses presented in the manuscript titled "Ribosome accumulation during early phase resistance training".

## Repository structure

Scripts are stored in the `./R` folder. `./data` contains all raw data with processed data (model outputs, summaries etc.) stored in `/data/derivedData`. `./resources` contains files used to produce outputs from R markdown files. All code creating figures (1-4) are located in the `./figures` folder together with resulting pdf-files.

The `manuscript.Rmd` file compiles the manuscript

### Western blot analyses

- western-compile.R: Compiles two western blot data sets, one allowing for comparisons of differences between groups over time (i.e. experimental vs. control and VAR vs. CONST). The second data set is a calibrated data set where participants in the experimental group are made comparable through calibration samples from each participant with the purpose of enabling analysis of the effect of UBF and S6 on total RNA.

- western-analysis.R: This script contains model formulations for brms fits. Comparisons of training effects (i.e. experimental vs. control) and volume effects (i.e. VAR vs. CONST).

### Total RNA analyses

- `total-rna-analysis.R`: Includes most analyses of total RNA. Data are compiled in `total-rna-compile.R`.

- `ubf-tot-rna-model.R`: UBF is believed to be a contributor to rDNA transcription by recruiting the SL1 transcription factor to the rDNA promoter. mTOR increases UBF phosphorylation (rapamycin sensitive) increasing its activity and c-Myc induces UBTF transcription increasing protein levels (rapamycin insensitive).

In this script we want to test if rRNA levels are associated with UBF levels. Both are increasing with time (western-analysis.R, total-rna-analysis.R) and therefore, the model controls for the effect of training (time).

This script also contain a model of muscle growth were we estimate the effect of increase in rRNA on muscle growth (muscle thickness change pre- to post training). Post training values are from both S12 and post 1w (de-training). rRNA transcription is estimated as the average linear increase over all sessions for each individual and leg.

### qPCR analyses

- `qpcr-compile.R` extracts cq-values and target-specific amplification efficiency from raw fluorescence data using the `qpcR` package.

- `qpcr-analysis-bayes2.R` models qPCR data in a Bayesian framework. Matz suggested using a poisson mixed-effects model to model qPCR-data without the use of reference genes as random effects on technical duplicates could be used to control for technical variation (Matz el al. 2013). An implementation of these models is done in the mcmc.qpcr package. Although this is a nice implementation we would like to have more control in the fitting process and downstream comparisons. To this end we use the `MCMCglmm` package directly. Two models are used, in "Model 1" estimates are retrieved from an "un-adjusted model" modeling transcript counts per total RNA (as this is the biological input). To account for technical variation in "Model 1", random effects are modeled for technical replicates as suggested by Matz et al. "Model 2" uses the external reference gene as a normalization factor modeling counts as a rate per tissue weight. The normalization factor is entered as a fixed effect with strong priors with the interpretation being counts per normalization factor. See discussion on offset in MCMCglmm: https://hannahdugdale.wordpress.com/2016/03/16/why-you-shouldnt-use-the-offset-function-in-mcmcglmm/ and https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q3/016817.html. 
 
### Training data, muscle mass and muscle strength

- `us-freq-bayes.R` fits models using ultrasound measurements. 

- `strength-bayes.R` fits models using strength measurements.

- `training-data.R` analyzes training data.

### Other R-files

`figure-source.R` contains color scales and other settings for figures. 

`libs.R` contains a list of all packages used in most scripts. Some scripts do not source this file to avoid to many packages being loaded. 



