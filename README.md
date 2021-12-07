# Ribosome accumulation during early phase resistance training

Daniel Hammarström<sup>1,2,£</sup>&ast;,Sjur J. Øfsteng<sup>2,£</sup>, Nicolai B. Jacobsen<sup>1</sup>, Krister B. Flobergseter<sup>1</sup>,  Bent R. Rønnestad<sup>1</sup>, Stian Ellefsen<sup>1,3</sup>


<sup>1</sup> Section for Health and Exercise Physiology, Department of Public Health and Sport Sciences, Inland Norway University of Applied Sciences, Lillehammer, Norway. 
<sup>2</sup> Swedish School of Sport and Health Sciences, Stockholm, Sweden.
<sup>3</sup> Innlandet Hospital Trust, Lillehammer, Norway. </div>
<sup>£</sup> These authors contributed equally to this work

----------------------------------------

This repository contains data, code and additional files needed to reproduce the analyses presented in the manuscript titeled "Ribosome accumulation during early phase resistance training". 

## Repository structure


### Western blot analyses

- western-compile.R: Compiles two western blot data sets, one allowing for comparisons of differences between groups over time (i.e. experimental vs. control and VAR vs. CONST). The second data set is a calibrated data set where participants in the experimental group are made comparable through calibration samples from each participant with the purpose of enabling analysis of the effect of UBF and S6 on total RNA.

- western-analysis.R: This script contains model formulations for brms fits. Comparisons of training effects (i.e. experimental vs. control) and volume effects (i.e. VAR vs. CONST).

### Total RNA analyses

- ubf-tot-rna-model.R: UBF is believed to be a contributor to rDNA transcription by recruiting the SL1 transcription factor to the rDNA promoter. mTOR increases UBF phosphorylation (rapamycin sensitive) increasing its activity and c-Myc induces UBTF transcription increasing protein levels (rapamycin insensitive).

In this script we want to test if rRNA levels are associated with UBF levels. Both are increasing with time (western-analysis.R, total-rna-analysis.R) and therefore, the model controls for the effect of training (time).

This script also contain a model of muscle growth were we estimate the effect of increase in rRNA on muscle growth (muscle thickness change pre- to post training). Post training values are from both S12 and post 1w (de-training). rRNA transcription is estimated as the average linear increase over all sessions for each individual and leg.

### Figures

All code creating figures (1-4) are located in the `./figures` folder.


