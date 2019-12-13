
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
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--12--13-yellowgreen.svg)](/commits/master)
[![](https://badges.ropensci.org/281_status.svg)](https://github.com/ropensci/software-review/issues/281)

**CytoExploreR** is comprehensive collection of interactive exploratory
cytometry analysis tools designed under a unified framework.
**CytoExploreR** has been specifically designed to integrate all
existing cytometry analysis techniques (e.g. manual gating, automated
gating and dimension reduction) in a format that makes these tools
freely accessible to users with minimal coding experience. If you are
new to **CytoExploreR** visit
<https://dillonhammill.github.io/CytoExploreR/> to get started.

## Install R and RStudio

**CytoExploreR** has a minimal requirement for R 3.5.0. If necessary,
newer versions of R can be installed by clicking on your operating
system in [this link](https://cran.r-project.org/bin/) and following the
installation instructions. For the best user experience it is
recommended that RStudio be installed as well. RStudio Desktop is free
to download and can installed from the
[RStudio](https://rstudio.com/products/rstudio/#rstudio-desktop)
website.

## Platform-Specific Requirements

After successfully installing R and RStudio, the following additional
platform-specific tools are required:

## Mac OS

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

## Windows OS

  - Install the appropriate
    [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R
    installation.
  - Follow these
    [instructions](https://github.com/RGLab/RProtoBufLib/blob/trunk/INSTALL)
    to download protobuf Windows binary and set appropriate environment
    variable. Ignore the section to build protobuf from source.
  - Follow these
    [instructions](https://github.com/RGLab/flowWorkspace/blob/trunk/INSTALL)
    to install and setup the additional C++ libraries required to
    successfully build flowWorkspace.
  - Restart your computer.

## CytoExploreR Dependencies

Now that all the setup is complete, let’s install all the necssary
dependencies of **CytoExploreR**.

From within RStudio, run the following in the console to install the
latest versions of [flowCore](https://github.com/RGLab/flowCore),
[flowWorkspace](https://github.com/RGLab/flowWorkspace) and
[openCyto](https://github.com/RGLab/openCyto) from Bioconductor.

``` r
# Bioconductor
install.packages("BiocManager")
# Install flowCore, flowWorkspace and openCyto
library(BiocManager)
install(c("flowCore", "flowWorkspace", "openCyto"))
```

Currently, **CytoExploreR** requires the development versions of these
RGLab packages hosted on GitHub:

``` r
# Install & load devtools
tryCatch(library(devtools), error = function(e){
  install.packages("devtools")
  library(devtools)
})
# Install flowCore, flowWorkspace & openCyto from GitHub
install_github("RGLab/flowCore", ref = "trunk")
install_github("RGLab/flowWorkspace", ref = "trunk")
install_github("RGLab/flowStats", ref = "trunk")
install_github("RGLab/openCyto", ref = "trunk")
```

## CytoExploreR

Now that all the dependencies are installed, let’s move on to installing
**CytoExploreR**.

To successfully install **CytoExploreR** users will first need to
install **CytoExploreRData** which contains example datasets that will
be used within **CytoExploreR** to demonstrate key features.

``` r
# CytoExploreRData 
devtools::install_github("DillonHammill/CytoExploreRData")
# CytoExploreR 
devtools::install_github("DillonHammill/CytoExploreR")
```

## Design

To ease the transition from GUI oriented software, **CytoExploreR** has
been designed to be a consistent and auto-complete friendly package for
cytometry data analysis. All exported functions from **CytoExploreR**
are prefixed with `cyto_`, followed by the name of the object you wish
to change (e.g. `cyto_gate_`) and finally the action that you would like
to perform (e.g. `cyto_gate_draw`). To see all available functions
simply start typing `cyto_` and you will be greeted with a complete list
of exported functions that can be selected from the auto-complete
dropdown
list:

<img src="man/figures/README-CytoExploreR-Design.gif" width="95%" style="display: block; margin: auto;" />

## Overview

Some of the key features of **CytoExploreR** are outlined below:

  - load and annotate samples using `cyto_setup`
  - user guided automatic compensation using `cyto_spillover_compute`
  - interactively modify spillover matrices using `cyto_spillover_edit`
  - compute spillover spreading matrices with
    `cyto_spillover_spread_compute`
  - visualise compensation in all channels using
    `cyto_plot_compensation`
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
  - exploratory visualisation of flowFrames, flowSets, GatingHierarchies
    and GatingSets using `cyto_plot`
  - visualisation of complete gating strategies with back-gating and/or
    gate tracking using `cyto_plot_gating_scheme`
  - visualisation of gating trees using `cyto_plot_gating_tree`
  - visualisation of marker expression profiles in all channels using
    `cyto_plot_profile`
  - visualisation of populations in all possible bivariate plots using
    `cyto_plot_explore`
  - produce dimension reduced maps (e.g. PCA, tSNE, UMAP and EmbedSOM)
    using `cyto_map`
  - save samples and analyses to file using `cyto_save`
  - export population level statistics tidyverse style using
    `cyto_stats_compute`

## Usage

The full details of how **CytoExploreR** works will be tackled
individually in the package vignettes, but a succinct usage outline is
described below:

**1. Create a new R project**

  - Select File -\> New Project
  - Select New Directory -\> New Project.
  - Set Directory Name to CytoExploreR-Demo.
  - Indicate where the project should be created and select Create
    Project.

**2. Save FCS files to new folders in the current working directory**

For demonstration purposes, let’s download the Compensation and
Activation FCS files shipped with CytoExploreRData:

``` r
# Load required packages
library(CytoExploreR)
library(CytoExploreRData)

# Download Compensation FCS files
cyto_save(Compensation,
          save_as = "Compensation-Samples")

# Download Activation FCS files
cyto_save(Activation,
          save_as = "Activation-Samples")
```

**3. Compensation of fluorescent spillover**

The basic steps to appropriate compensate for fluorescent spillover are
outlined below. Refer to the Compensation vignette for a more detailed
workflow.

  - Load Compensation-Samples into a GatingSet and prepare
    gatingTemplate

<!-- end list -->

``` r
gs <- cyto_setup("Compensation-Samples",
                 gatingTemplate = "Compensation-gatingTemplate.csv")
```

  - Gate compensation controls using `cyto_gate_draw`

<!-- end list -->

``` r
# Cells
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))

# Single Cells
cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A", "FSC-H"))
```

<img src="man/figures/README-Compensation-Gating.png" style="display: block; margin: auto;" />

  - Automated computation of spillover matrix using
    `cyto_spillover_compute`

<!-- end list -->

``` r
spill <- cyto_spillover_compute(gs,
                                parent = "Single Cells",
                                spillover = "Spillover-Matrix.csv")
```

  - Interactively edit computed spillover matrix using
    `cyto_spillover_edit`

<!-- end list -->

``` r
spill <- cyto_spillover_edit(gs,
                             parent = "Single Cells",
                             spillover = "Spillover-Matrix.csv")
```

<img src="man/figures/README-Compensation-Editing.png" style="display: block; margin: auto;" />

**4. Load and annotate Activation samples**

The basic steps to appropriately compensate, transform and manually gate
the Activation samples are outlined below. Refer to the Manual Gating
vignette for a more detailed workflow.

  - Load Activation-Samples into a GatingSet and prepare gatingTemplate

<!-- end list -->

``` r
gs <- cyto_setup("Activation-Samples",
                 gatingTemplate = "Activation-gatingTemplate.csv")
```

  - Apply compensation of fluorescent spillover using `cyto_compensate`

<!-- end list -->

``` r
gs <- cyto_compensate(gs, "Spillover-Matrix.csv")
```

  - Transform fluorescent channels using `cyto_transform`

<!-- end list -->

``` r
gs <- cyto_transform(gs, 
                     type = "logicle")
```

  - Manually gate samples to build gating scheme using `cyto_gate_draw`

<!-- end list -->

``` r
# Cells
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))

# Copy and paste above, modify parent and alias sequentially to add new populations
```

**5. Visualise gating scheme using `cyto_plot_gating_scheme`**

``` r
cyto_plot_gating_scheme(gs[[32]], back_gate = TRUE)
```

**6. Export population-level statistics using `cyto_stats_compute`**

``` r
cyto_stats_compute(gs,
                   alias = c("CD4 T Cells", "CD8 T Cells"),
                   channels = c("CD44", "CD69"),
                   stat = "median")
```

**7. Save GatingSet for future analyses**

``` r
cyto_save(gs,
          save_as = "Activation-GatingSet")
```

## News

There is a Changelog for the GitHub `master` branch which will reflect
any updates made to improve the stability, usability or plenitude of the
package. Users should refer to the
[Changelog](https://dillonhammill.github.io/CytoExploreR/news/index.html)
before installing new versions of the package.

## Credits

**CytoExploreR** would not be possible without the existing flow
cytometry infrastructure developed by the RGLab. **CytoExploreR**
started out as simple plugin for openCyto to facilitate gate drawing but
has evolved into a fully-fledged cytometry analysis package thanks to
the support and guidance of members of the RGLab. Please take the time
to check out their work on [GitHub](https://github.com/RGLab).

## Development

**CytoExploreR** is a maturing package which will continue to be
sculpted by the feedback and feature requests of users. The GitHub
`master` branch will always contain the most stable build of the
package. New features and updates will be made to a separate branch and
merged to the `master` branch when stable and tested. The
[Changelog](https://dillonhammill.github.io/CytoExploreR/news/index.html)
will reflect any changes made to the `master` branch.

## Getting help

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

## Code of conduct

Please note that the **CytoExploreR** project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to
this project, you agree to abide by its terms.
