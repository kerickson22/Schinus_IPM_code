# Rapid maturation of seedlings drives increased population growth rate of an invader with multiple introductions: insights from an LTRE analysis

# Kelley D. Erickson, Paul D. Pratt, Min B. Rayamajhi and Carol C. Horvitz

# Please contact Kelley D. Erickson for questions about the code or analysis (kerickson22@gmail.com)

# Abstract 

Multiple introductions are hypothesized to facilitate the success of invasive plant species, as these introductions can result in novel genotypes through intraspecific hybridization and potentially increasing the ability to adapt to the novel environment. Through this combination of genes and the environment, how an invader performs may vary across its invaded range. The ability to evaluate the success of any management action such as the introduction of a biocontrol agent could therefore be spatially dependent and requires an understanding of the population dynamics of the invasive species prior to management action. Brazilian Pepper (*Schinus terebinthifolia* Raddi), a shrub that has invaded the global subtropics, was introduced in Florida in two separate introductions (an Eastern biotype and a Western biotype) and intraspecific hybridization has resulted in a hybrid biotype. We constructed integral projection models where the probabilities of survival, growth and reproduction were functions of two size variables (diameter and height). To decompose the effects of the different biotypes on the population dynamics of Brazilian Pepper in the absence of any biocontrol we performed a Life Table Response Experiment analysis and found that the higher population growth rate of the Western population was driven by the comparatively rapid maturation of the Western seedlings into reproductive adults. 

# Code 
* 1demography_models.R : This script performs the statistical analyses to model the vital rates of survival, growth, reproduction, fecundity and graduation from the D1 size domain to the D2 size domain. 
* 2IPM_Diameter_and_Height.R : This script constructs the integral projection models from the demographic statistical models. 
* 3LTRE_analysis.R: This script performs a Life Table Response Experiment analysis using the integral projection models constructed from the previous script. 
* 4MSFigures.R : This script creates the publication-ready figures for the manuscript. 
