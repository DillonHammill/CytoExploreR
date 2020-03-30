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