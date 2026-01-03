# Hadlock-Growth-Chart-Statistical-Validation-Analysis
THE HADLOCK REFERENCE CHART FOR FETAL WEIGHT: TIME TO RECONSIDER ITS CLINICAL USE?

**Authors**: Theophilus Adu-Bredu, Clement Twumasi, Sam Mathewlynn, Christos Ioannou, Lawrence Impey

# Internal Consistency of the Hadlock Fetal Growth Reference

# R code repository

**Overview:**

This repository contains the R code (dubbed `Main R codes_HadlockChart_Validation.r` and `Empirical-Evaluation-Hadlock-Centiles.r`) used to assess the internal consistency between the published Hadlock fetal growth reference charts and the corresponding equation-derived centiles, and to evaluate the clinical implications for classification of small for gestational age (SGA) fetuses and perinatal outcomes.

The analyses underpin the findings reported in the associated manuscript, including discrepancies in centile classification (particularly the 3rd and 10th centiles), agreement between chart-based and equation-based approaches, and associations with adverse perinatal outcomes.

# Objectives Addressed by the Code

1. Evaluate the agreement between the Hadlock reference charts and the Hadlock equation-derived centiles

2.  Quantify misclassification of SGA at clinically relevant thresholds (3rd and 10th centiles)

3. Apply functional data analytic methods to assess internal consistency across gestation


# Key Analytical Components

1. Hadlock Equation Implementation

Reproduces the original Hadlock fetal weight equations

Generates equation-derived centiles across gestational age

Ensures numerical consistency with published parameters

2. Chart-Based vs Equation-Based Centiles

* Digitised interpolation of published Hadlock reference charts

* Direct comparison of:

  - 3rd centile

  - 10th centile

* Identification of centile drift between methods

3. Functional Data Analysis

* Treats centile curves as continuous functions

* Quantifies structural inconsistency between the chart and the equation

* Highlights gestational-age–specific divergence


# Description of the R scripts:

1. `Main R codes_HadlockChart_Validation.r`: This file contains all the R codes used to support the numerical analysis (including the functional data analysis and other statistical analysis) on the internal consistency between the standard Hadlock growth chart and the corrected growth chart (derived in the current study) under Hadlock's optimal growth model.
2. `Empirical-Evaluation-Hadlock-Centiles.r`: This file contains the R codes used to support additional empirical statistical evaluation of the Hadlock growth chart discrepancy at some pre-selected week gestational ages and its corresponding estimated fetal weight (EFW) centiles based on an observed/empirical clinical (estimated fetal weight) dataset. 

**NB: The clinical dataset used for the statistical analyses is available from the corresponding author upon reasonable request.**

