# SPACE
Spatial Patterning Analysis of Cellular Ensembles (SPACE), an R Package

## Overview
SPACE is an R package for spatial analysis of multiplex images of biological tissues, although spatial data from any context can be used. For many ensembles (i.e. combinations) of spatial elements (e.g. cell types, biomolecules, etc.), SPACE quantifies non-random patterning in a single specimen (called "cisMI") or divergent patterning across groups of specimens (called "transMI"). Thanks to the information-theoretic underpinnings of SPACE, these patterns can be arbitrarily complex: simple co-assortment and mutual exclusion are detected, along with quantitative gradients, structural orientations, context-dependent interactions of 3+ elements, etc. Further tools allow in-depth characterization and visualization of specific ensembles. 

## Contents and Demo
In addition to the source R code, this package also includes:
- documentation detailing the purpose, inputs, parameters, and outputs of each function
- a tutorial describing a typical SPACE analysis workflow
- sample data to follow along with the tutorial

To get started, follow the tutorial, which includes all commands, suggested parameters and arguments, notes on expected runtimes, and expected outputs.

## Current Status
A manuscript that details and demonstrates SPACE is posted on bioRxiv (https://www.biorxiv.org/content/10.1101/2023.12.08.570837v2) and is also currently under review for publication.

## Installation
SPACE can be installed directly from Github:
- Install the R package "devtools" if you don't already have it.
- Run the command devtools::install_github("eschrom/SPACE")
- Restart your R session

Installation should only take a few minutes, but perhaps more if many of the dependencies must also be installed or updated.

After installing SPACE, to access the tutorial and its associated materials, locate the SPACE folder and open the doc subfolder. "SPACE_Tutorial.html" is the tutorial itself, and the other files in the doc subfolder are the supporting data.

## Computing Requirements and Runtimes
SPACE was developed on a PC laptop with a 11th Gen Intel(R) Core(TM) i9-11900H (2.50 GHz) processor and 16.0 GB RAM, running Windows 10 64-bit operating system, and using R version 4.3.1 "Beagle Scouts." On this system, most SPACE functions operate in seconds to minutes or faster. However, the 'census_image', 'measure_cisMI', and 'measure_transMI' functions can require several hours, depending on the size of the data and number of calculations requested. The minimum computing requirements are smaller than listed here, but compute times will scale with the size of the data and your machine's capacity.
