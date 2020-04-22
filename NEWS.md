# CytoExploreR 1.0.7

* Fix quadrant gating bug in `cyto_gate_draw()`.
* Switch to using `rsvd::rpca()` for faster PCA computation in `cyto_map()`.
* Add versioned and latest tags to CytoExploreR Docker images for reproducibility. Docker images are always built against the most recent RGLab packages at the time. Versioned Docker images allow users to go back in time should they need to.
* Update FIt-SNE to version `1.2.1`. This change requires updated installation of FIt-SNE and now requires the `rsvd` package. FIt-SNE now uses PCA initialization by default.
* Configure Docker images with FIt-SNE and search for `fast_tsne` location automatically. Hide this directory to prevent accidental deletion.
* Add FIt-SNE configuration instructions on GitHub [Issue #29](https://github.com/DillonHammill/CytoExploreR/issues/29). Configuration instructions for mac OS are still a work in progress.
* Make improvements to flowAI package to improve compatibility with CytoExploreR. FlowAI can be used when loading samples in `cyto_setup()` by setting `clean = TRUE` or manually applied downstream using `cyto_clean()`.
* Add back Appveyor CI to build on windows.

# CytoExploreR 1.0.6

* Update custom plotting functions (`cyto_plot_save`, `cyto_plot_layout`, `cyto_plot_custom` and `cyto_plot_complete`) to allow complex layouts. Layouts must now be supplied to these functions either as a vector or matrix.
* `cyto_plot_complete` no longer closes the RStudio or X11 graphics devices.
* Add new `margins` argument to `cyto_plot` to allow customization of plot margins.
* Fix `overlay = "children"` when plotting GatingSet objects with `cyto_plot`.

# CytoExploreR 1.0.5

* Add ability to control plotting order in `cyto_plot` when using the `group_by` argument.
* Add new `data_editor` function to bring excel-like features to all context menus and data editors.
* Add viewer argument to `cyto_spillover_edit` to allow it to be launched in the RStudio viewer pane (if desired).

# CytoExploreR 1.0.4

* Remove requirement for GUI device in `cyto_gate_draw` and `cyto_gate_edit`. Users can now draw gates directly on the RStudio graphics device.
* Add support for changing gate point and lines colours within `cyto_gate_draw` and `cyto_gate_edit`.
* improve appearance of 1-D interval gates in `cyto_gate_draw`.
* Add dependency for openCyto >= 1.25.2 to support gating of empty populations.
* Add Dockerfile and Makefile to build **CytoExploreR** docker images for reproducible data analysis. Docker images can be pulled from [docker hub](https://hub.docker.com/) at `dhammill/cytoexplorer`.

# CytoExploreR 1.0.3

* Update installation instructions on website following in changes in RGLab default branches.
* Fix sample grouping when no experiment details are added to samples.
* Resort to machine axes limits when no events are present for plotting. 
* Improve handling of multiple compensation controls for the same channel in `cyto_spillover_compute()`. Now the control with the best signal will be used for the calculation.
* Add support to `cyto_map()` for FIt-SNE dimensionality reduction. Update `cyto_map()` documentation to include references for each of the supported dimensionality reduction algorithms.
* Fix gating of empty populations to prevent dummy gates being returned by openCyto.
* Frequency statistics for empty populations are now displayed as zero instead of NaN.
* Improve handling of arguments within cyto_plot for better user experience.
* Fix errors for cyto_spillover_compute and cyto_spillover_spread_compute on mac OS.
* `cyto_setup()` now accepts filenames through markers and details arguments.

# CytoExploreR 1.0.0

* Added a `NEWS.md` file to track changes to the package.