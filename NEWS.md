# CytoExploreR 1.0.9

* Add new `cyto_cbind()` function to handle addition of new parameters to cytometry objects (retaining support for flowFrames/flowSets).
* `cyto_apply()` has been re-written to accept different data inputs (e.g. matrix, cytoframe, cytoset, column or row) and the output is more intuitively formatted.
* `cyto_calibrate()` now uses quantile calibration by default and the upper quantile limits has been reduced from 0.99 to 0.95.
* Barcoding of GatingSets is now supported by `cyto_barcode()`.
* `cyto_beads_sample()` is now defunct in favour of the more flexible `cyto_sample_to_node()` function.
* Add new `cyto_sort_by()` function to sort samples based on experimental variables.
* Add new `cyto_groups()` function to return the names or details experimental groups.
* `cyto_extract()` is now defunct in favour of the new `cyto_data_extract()` which has for extracting multiple populations. New arguments `format`, `split`, `markers`, `trans` and `inverse` provide flexibility in how the extracted population should be formatted.
* `cyto_names()` has been updated to extract sample information at the level of cytosets only, attempts to extract this information from flowFrames/cytoframes will generate an error.
* `cyto_details()` now has additional `convert`, `drop` and `factor` arguments to allow easy formatting of the returned data.frame. Support for extracting `cyto_details()` from flowFrame/cytoframe objects has been removed.
* `cyto_transform()` now appropriately handles transforming subsetted data where some parameters defined in the transformerList may be missing in the data.
* `cyto_details()` called on a cytoframe now returns a data.frame to match the behaviour of cytosets, GatingHierarchies and GatingSets.
* Add new `cyto_class()` to get or check the class of cytometry objects.
* Remove directories from file paths in `cyto_load()`. Re-factor code to ensure loading of saved GatingSets.

# CytoExploreR 1.0.8

* Add `select` argument to `cyto_stats_compute()` to allow computation of statistics for a subset of samples based on experimental variables.
* Add new `cyto_gate_bool()` function to add boolean gates to the GatingSet and gatingTemplate.
* Improve `cyto_channels_restrict()` to drop unused channels in cytometry objects. Add an `exclude` argument to forcibly remove FSC, SSC or Time parameters. The `restrict` argument in `cyto_setup()` now accepts a vector of channels to be passed internally to the `exclude` argument of `cyto_channels_restrict()`.
* Add support for forward gating by colouring points based on the expression of a particular marker in `cyto_plot()`, simply pass the name of the marker to the `point_col` argument. 
* Update default density colour scale for points to include a darker blue at the lower end of the scale. This provides better resolution at the lower end of the scale and balances the colour gradient.
* Add channels argument to `cyto_extract()` to restrict data to certain channels when it is extracted from a GatingHierarchy or GatingSet.
* Add new function `cyto_names_parse()` to split file names into experiment variables in `cyto_details()`. Files should be named with each variable separated by a delimiter. For example, file named 160520-Exp1-StimA.fcs can be separated into variables date, experiment and treatment. `parse_names` argument has been added to `cyto_setup()` to allow name parsing when files are loaded in.
* Custom functions can now be supplied by name to the `type` argument in `cyto_map()`. This means that you can use any available dimensionality reduction technique or create your own one to map your data.
* Add `cyto_gatingTemplate_generate()` to obtain a CytoExploreR-ready gatingTemplate from any GatingHierarchy or GatingSet.
* Improve speed of plotting using `cyto_plot()` by increasing computation speed of axes ranges.
* Add support for `scattermore` to further improve `cyto_plot()` plotting speed. To use this feature, users will need to manually install `scattermore` and set `options("cyto_plot_fast" = TRUE)`. Setting `options("cyto_plot_fast" = FALSE)` will return to using conventional plotting. It is recommended to use conventional plotting when saving images with `cyto_plot_save()` to obtain high resolution images.
* Fix quadrant gate bug in `cyto_gate_edit()` to allow editing of quadrant gates.
* Fix computation of statistics for stacked density distributions in `cyto_plot()`.
* Fix `cyto_marker_edit()` to accept custom file name for saving and add support for marker removal.
* Hide messages when opening data editors and only open data editors in interactive mode.
* Improve visualization of quadrant gates to distinguish them from traditional rectangleGates. Quadrant gates can now be individually visualized (i.e. users can indicate which quadrants to display).
* Added new function `cyto_nodes_convert()` to anchor ambiguous nodes to a known parental node within `cyto_plot()`. 

# CytoExploreR 1.0.7

* Fix quadrant gating bug in `cyto_gate_draw()`.
* `cyto_map()` can now appropriately handle empty flowFrames/cytoframes.
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