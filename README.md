
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CytoExploreR <img src="man/figures/logo.png" align="right" alt="" width="130"/>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Travis build
status](https://travis-ci.org/DillonHammill/CytoExploreR.svg?branch=master)](https://travis-ci.org/DillonHammill/CytoExploreR)
[![Coverage
status](https://codecov.io/gh/DillonHammill/CytoExploreR/branch/master/graph/badge.svg)](https://codecov.io/github/DillonHammill/CytoExploreR?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--11--09-yellowgreen.svg)](/commits/master)
[![](https://badges.ropensci.org/281_status.svg)](https://github.com/ropensci/software-review/issues/281)

**CytoExploreR** is comprehensive collection of exploratory cytometry
analysis tools under a unified interactive framework. If you are new to
**CytoExploreR** visit <https://dillonhammill.github.io/CytoExploreR/>
to get started.

# Installation

## R

**CytoExploreR** has a minimal requirement for R 3.5.0. If necessary,
newer versions of R can be installed by clicking on your operating
system in [this link](https://cran.r-project.org/bin/) and following the
installation instructions.

## RStudio

For the best user experience it is recommended that RStudio be installed
as well. RStudio Desktop is free to download and can installed from the
[RStudio
website](https://rstudio.com/products/rstudio/#rstudio-desktop).

After successfully installing R and RStudio, the following additional
platform-specific tools are required:

## Mac

  - Install Xcode developer tools from the [App
    Store](https://apps.apple.com/au/app/xcode/id497799835?mt=12).
    Restart your computer.
  - Install command line tools by opening the terminal and running
    xcode-select –install
  - Install macOS R toolchain by installing clang7 and gfortran
    [here](https://cran.r-project.org/bin/macosx/tools/)
  - Check if XQuartz is listed in your installed Applications, it may be
    hiding in the utilities folder. If XQuartz is missing on your
    computer, it can be installed from the
    [XQuartz](https://www.xquartz.org/) website.
  - Restart your computer so all these changes will take effect.

## Windows

  - Install the appropriate
    [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R
    installation.
  - Follow these
    [instructions](https://github.com/RGLab/flowWorkspace/blob/trunk/INSTALL)
    to install and setup the additional C++ libraries required to
    successfully build flowWorkspace.
  - Restart your computer.

Now that all the setup is complete, let’s install all the necssary
dependencies of **CytoExploreR**:

## flowCore, flowWorkspace & openCyto

From within RStudio, run the following in the console to install the
latest versions of[flowCore](https://github.com/RGLab/flowCore),
[flowWorkspace](https://github.com/RGLab/flowWorkspace) and
[openCyto](https://github.com/RGLab/openCyto) from Bioconductor.

``` r
# Bioconductor
install.packages("BiocManager")
# Install flowCore, flowWorkspace and openCyto
library(BiocManager)
install(c("flowCore", "flowWorkspace", "openCyto"))
```

Now that all the dependencies are installed, let’s move on to installing
**CytoExploreR**:

## CytoExploreR

To successfully install **CytoExploreR** users will first need to
install **CytoExploreRData** which contains example datasets that will
be used within **CytoExploreR** to demonstrate key features.

``` r
# CytoExploreRData development version on GitHub
devtools::install_github("DillonHammill/CytoExploreRData")
# CytoExploreR development version on GitHub
devtools::install_github("DillonHammill/CytoExploreR")
```

# Overview

**CytoExploreR** provides an interactive interface for analysis of flow
cytometry data. Some key features include:

  - user guided automatic compensation using `cyto_spillover_compute`
  - interactively modify spillover matrices using `cyto_spillover_edit`
  - compute spillover spreading matrices with
    `cyto_spillover_spread_compute`
  - visualise compensation in all channels using
    `cyto_plot_compensation`
  - easily associate experimental details with each file using
    `cyto_annotate`
  - customisable data transformations using `cyto_transform` which
    includes support for log, arcsinh, logicle and biexponential
    transformations
  - manual gate drawing using `cyto_gate_draw`
  - ability to edit drawn gates using `cyto_gate_edit`
  - remove gates using `cyto_gate_remove`
  - rename gates using `cyto_gate_rename`
  - gate saving directly to an openCyto `gatingTemplate` for future use
  - support for using both manual and automated gating approaches
    through linking to `openCyto`
  - visualisation of flowFrames, flowSets, GatingHierarchies and
    GatingSets using `cyto_plot`
  - visualisation of complete gating strategies with back-gating and/or
    gate tracking using `cyto_plot_gating_scheme`
  - visualisation of gating trees using `cyto_plot_gating_tree`
  - visualisation of marker expression profiles in all channels using
    `cyto_plot_profile`
  - visualisation of populations in all possible bivariate plots using
    `cyto_plot_explore`
  - export population level statistics tidyverse style using
    `cyto_stats_compute`

# Usage

# News

There is a Changelog for the GitHub `master` branch which will reflect
any updates made to improve the stability, usability or plenitude of the
package. Users should refer to the
[Changelog](https://dillonhammill.github.io/CytoExploreR/news/index.html)
before installing new versions of the package.

# Credits

**CytoExploreR** would not be possible without the existing flow
cytometry infrastructure developed by the RGLab. **CytoExploreR**
started out as simple plugin for openCyto to facilitate gate drawing but
has evolved into a fully-fledged flow cytometry analysis package thanks
to the support and guidance of members of the RGLab. Please take the
time to check out their work on [GitHub](https://github.com/RGLab).

# Development

**CytoExploreR** is a maturing package which will continue to be
sculpted by the feedback and feature requests of users. The GitHub
`master` branch will always contain the most stable build of the
package. New features and updates will be made to a separate branch and
merged to the `master` branch when stable and tested. The
[Changelog](https://dillonhammill.github.io/CytoExploreR/news/index.html)
will reflect any changes made to the `master` branch.

# Getting help

The [Get
Started](https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR.html)
and
[Reference](https://dillonhammill.github.io/CytoExploreR/reference/index.html)
sections on the **CytoExploreR** website are your first port of call if
you require any help. For more detailed workflows refer the **Articles**
tab. If you encounter any issues with the functioning of the package
refer to these
[issues](https://github.com/DillonHammill/CytoExploreR/issues) to see if
the problem has been identified and resolved. Feel free to post new
issues on the GitHub page if they have not already been addressed.

# Code of conduct

Please note that the **CytoExploreR** project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to
this project, you agree to abide by its terms.
