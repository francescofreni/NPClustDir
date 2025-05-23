# Nonparametric Clustering for Directional Data

This repository contains the code for the conference paper *Efficient Disentangling of &#947;-Ray Sources from Diffuse Background in the Sky Map*, presented at the 14th Scientific Meeting of the Classification and Data Analysis Group (CLADAG) of the Italian Statistical Society (SIS) in Salerno, Italy, on September 11-13, 2023. You can find the short paper [here](https://it.pearson.com/content/dam/region-core/italy/pearson-italy/pdf/Docenti/Universit%C3%A0/CLADAG-2023.pdf).

## Contributors
- Francesco Freni
- Giovanna Menardi

## Abstract
Searching for as yet undetected &#947;-ray sources is a major target of the Fermi LAT Collaboration. We address the problem by clustering the directions of the high-energy photon emissions detected by the telescope onboard the Fermi spacecraft. Putative sources are identified as the excess mass of disconnected high density regions on a sphere mesh, which allows for their joint discrimination from the diffuse &#947;-ray background spreading over the entire area.
Density is estimated nonparametrically via binned directional kernel methods. The identification is accomplished by breaking the problem into independent subregions of the sphere separated by empty bins, thus leading to a remarkable gain in efficiency.

## Getting Started
To get started with the project, clone the repository by running the following command:
```bash
git clone https://github.com/francescofreni/fermi-clust.git
cd fermi-clust
```

## Usage
The scripts are organized as follows:

- **`R/`**: Contains all the R scripts, including function definitions and a simple example.
- **`R/healpix/`**: Includes adaptations of the functions to work with the HEALPix pixelization.

