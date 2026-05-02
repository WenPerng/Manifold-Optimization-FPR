# Manifold Optimization on the Magnitude Torus for Fourier Phase Retrieval

> This repository contains the code to reproduce the simulations present in our work "*Manifold Optimization on the Magnitude Torus for Fourier Phase Retrieval*," which includes the code used in both the manuscript and the supplemental materials.

## Descriptions of Codes
Here we separate the codes by their purposes, so the user can read the descriptions and know which one they want to use and download.

The algorithms we provided (replaced with `XXX` below) includes:
- the Gerchberg--Saxton algorithm (`GS`), 
- the projected gradient descent algorithm (`PGD`),
- the continuous hybrid input-output algorithm (`CHIO`),
- the dualPGD algorithm (`dualPGD`), and
- the TorusGD algorithm (`TorusGD`).

### Success Rate
The codes with `SuccessRate_XXX` in their filenames are used to generate the success rate seen in **Table 1** of the paper, which outputs the success rate of the said algorithm under 1000 iterations.

### Algorithm
The codes with `Algorithm_XXX` in their filenames simply runs the specified algorithm.

### Data Set
The file `img_BadApple.mat` is a saved image variable that loads the test image **Fig. 3(d)**.