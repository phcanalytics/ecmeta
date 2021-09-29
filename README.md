# Meta-analysis for external control studies
## Overview
The `ecmeta` `R` package provides tools to implement a meta-analytic framework that adjusts external control studies for additional bias and variability due to their non-randomized design. The meta-analysis is used to compare internal trial controls and external controls and is performed on a set of historical references studies. The parameters of the meta-analytic model are, in turn, used to to adjust estimate of a new study that compares an experimental treatment with an external control arm.

A paper describing the methodology and evauating its performance is available at the GitHub repository [here](https://github.roche.com/incertid/ecmeta-manuscript). For an example analysis, see the [user guide](https://pages.github.roche.com/RWDScodeshare/ecmeta/articles/guide.html).

## Installation
### R package
The `R` package can be installed with:

```{r}
devtools::install_git("https://github.roche.com/RWDScodeshare/ecmeta.git")
```

### JAGS
The meta-analytic models can be estimated using Bayesian techniques with `JAGS`. If you would like to use `JAGS`, you must first install it. On Linux:

```
sudo apt update
sudo apt-get install jags
```

On OS X or Windows, you can install from the [website](https://mcmc-jags.sourceforge.io/).

Afterwards, you will also need to install the `rjags` package, which provides an `R` based interface:

```
install.packages("rjags")
```
