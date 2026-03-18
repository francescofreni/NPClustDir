# Nonparametric Clustering for Directional Data

This repository contains the code for the research paper [*Efficient Disentangling of &#947;-Ray Sources from Diffuse Background in the Sky Map*](https://www.sciencedirect.com/science/article/pii/S2590197426000200).

## Contributors
- Francesco Freni
- Giovanna Menardi

## Abstract
Searching for previously unknown &#947;-ray sources is a key objective of the Fermi LAT Collaboration. We address this challenge by clustering the directions of high-energy photon emissions detected by the Fermi spacecraft’s Large Area Telescope (LAT). Candidate sources are detected by analyzing the excess mass within discrete, high-density regions, allowing us to discriminate them from the diffuse &#947;-ray background that pervades the entire sky. Density estimation is performed nonparametrically using binned directional kernel methods applied to a sphere mesh. Source detection is facilitated by partitioning the problem into separate subregions of the sphere, delimited by empty bins, which results in a substantial gain in computational efficiency.

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

