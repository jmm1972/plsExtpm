# plsExtpm
R code package for Partial Least Squares Path Modelling allowing estimation of nonlinear structural relationships

Partial Least Squares Path Modelling is an iterative method used to estimate Structural Equation Models (SEM), a widely used analytical tool for assessing causal relationships between latent variables.

A SEM consists of two models: a measurement model (outer model) and a structural model (inner model). While the former groups manifest variables (observable variables or indicators) into corresponding latent variables (or constructs), the latter comprises a set of regression equations assessing the effects of explanatory latent variables on the dependent latent variable. Over the past 20 years, SEM, also referred to as a second-generation multivariate technique, has gained popularity.

However, PLS-PM struggles to address the structural nonlinear relationships. To address this limitation, a new PLS-PM inner weighting scheme, smooth weighting, is proposed as an additional option to the traditional centroid, factor, and path weighting schemes. This repository stores the code I developed that will be continuously improved.

It also contains the article XXX and the description of the code used to produce the article's results: (1) a relationship of interest in marketing, the relationship between customer satisfaction and customer loyalty and (2) a simulated dataset is used to assess the ability of the algorithm to approximate the underlying (unknown) nonlinear structural relationships. The results show that the proposed scheme can recover several nonlinear functional forms, outperforming existing inner weighting schemes for commonly used sample sizes (larger than 75 units), regardless of the level of error contamination in the observed manifest variables.

if you want to use it, feel free! But don't forget to cite the work:

To cite package ‘plsExtpm’ in any publication, use:

  Mendes, J. M. and Coelho, P. S. (2024). plsExtpm: Extended Partial Least Squares Path Modeling (PLS-PM). An R
  unofficial package, https://github.com/jmm1972/plsExtpm.

A BibTeX entry for LaTeX users is

  @article{,
    title = {plspm: Partial Least Squares Path Modeling (PLS-PM)},
    author = {Jorge M. Mendes and Pedro S. Coelho},
    year = {2024},
    journal = {to be announced},
    url = {https://github.com/jmm1972/plsExtpm},
  }

The official introduction of the package will be done through a scientific article to be submitted soon.
