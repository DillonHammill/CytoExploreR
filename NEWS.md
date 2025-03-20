# CytoExploreR 2.0.13 (pre-release)

* The behavior of `cyto_merge_by()` when `merge_by = NA` has changed from collapsing all samples to instead split samples individually. This is because splitting by `name` may not always work in cases where multiple samples share the same file names. This change requires updates to openCyto which should be re-installed when updating to the new version of CytoExploreR.

# CytoExploreR 2.0.12

* Group progress messages for `cyto_gate_draw()`
* Improved handling of `ylim` for mixed 1D/2D `cyto_plot()` wrappers.
* Consistent titles for all panels in `cyto_plot_explore()`
* Default histogram colour has changed to match that of points - the bottom colour of the density colour scale.
* Histogram and point colours are now in sync.
* `hist_fill` and `hist_fill_alpha` inherit from `point_col` and `point_col_alpha` to ensure consistency in multi-plot layouts. 

# CytoExploreR 2.0.11

* Bug fix to ensure overalys are binned correctly by default when `point_shape = "hex"`.

# CytoExploreR 2.0.10

* Added support for non-alphabet characters in `pData` grouping variables for `cyto_gate_draw()` openCyto plugin.

# CytoExploreR 2.0.9

* Added a more informative error message for use of overlays within `cyto_plot()`.

# CytoExploreR 2.0.8

* Fix handling of `multiRangeGate` class in `cyto_gate_edit()` to allow 1-D editing of `multiRangeGates`.

# CytoExploreR 2.0.7

* Bug fixes to handling of `page_fill` and `page_fill_alpha` within the `cyto_plot()` family of functions.
* Add `close` button to `cyto_spillove_edit()` which saves the spillover matrix and then kills the application.

# CytoExploreR 2.0.6

* `cyto_plot_calibrate()` can now use the full instrument ranges for calibration through `type = "instrument"`.

# CytoExploreR 2.0.5

* `cyto_import()` has gained a logical `nodes` argument that allows for bypassing of the unique nodes test and renaming (it can be temperamental for complex node names).
* `cyto_plot()` gain additional customization arguments including `axes_ticks_line_width`, `axes_ticks_line_col`, `axes_ticks_line_col_alpha`, `key_border_line_width`, `key_border_line_col`, `key_border_line_col_alpha`, `key_ticks_line_width`, `key_ticks_line_col`, `key_ticks_line_col_alpha`, `page_fill` and `page_fill_alpha`.
* `cyto_plot()` now has native support for hex binning through `point_col = "hex"` and binning can be controlled through `point_bins`. `point_bins` also influences the 2D kernel density estimates for points to allow for finer control over the resolution of the binned kernel density estimate. Hex bins are computed using the new `cyto_stat_hex()`.
* New custom themes have been added to match the theme of Ozette's platform including `ozette_dark_theme()` and `ozette_light_theme()`.

# CytoExploreR 2.0.4

* Gates can now be skipped upon review in `cyto_gate_edit()` by hitting `Finish` before selecting any points.

# CytoExploreR 2.0.3

* CytoExploreR can now read in spillover matrices from `.mtx` files using the new internal `read_from_mtx()` function.
* CytoExplreR can now be run without using a gatingTemplate by setting `cyto_gatingTemplate_active(FALSE)` at the beginning of the session.
* `cyto_plot_map()` now supports plotting gates through the `alias` argument to maintain consistency with `cyto_plot()`.
* Improvements to consistency of return values of `cyto_channels()`, `cyto_markers()` and associated functions.

# CytoExploreR 2.0.2

* Major improvements to `cyto_gatingTemplate_generate()` to allow for finer control over how gatingTemplates should be generated from GatingHierarchy or GatingSet objects. 
* `cyto_groups()` now reports the number of samples per group when `details = FALSE`.
* Changed default value of `multiple` to TRUE in `cyto_plot_save()`.

# CytoExploreR 2.0.1

* Added new `multirange` gate type to `cyto_gate_draw()` to allow for Time channel gating.
* `cyto_plot()` now supports visualisation of `multiRangeGate` objects.

# CytoExploreR 2.0.0

* `CytoExploreR` now uses the `HeatmapR` package to display heatmaps in `cyto_spillover_compute()` and `cyto_spillover_spread_compute()`.
* Added new `cyto_gate_clust()` gating function to apply ANY clustering algorithm and store populations as gates in the GatingHierarchy or GatingSet.
* Added new `cyto_gate_clean()` gating function to apply ANY anomaly detection algorithm and store high quality events as a gated population.
* Added new `cyto_gate_sample()` gating function to add sampled nodes to a GatingHierarchy or GatingSet objects. `cyto_gate_sample()` performs groupwise proportionate sampling by making calls to `cyto_sample_n()` on each experimental group. 
* `cyto_transformers_define()` has been updated to automatically compute parameter estimates for instruments with larger dynamic ranges, allowing for appropriate transformation of data from newer instruments.
* `cyto_channel_match()` can now automatically detect the channel associated with each compensation control by matching marker and channel combinations to the file names. Samples where channels cannot be matched by file name will instead be matched by the intensities in the fluorescent channels. `cyto_channel_match()` also gains the ability to automatically detect which compensation controls are beads or cells, an important distinction for `cyto_spillover_compute()`.
* Added new `label_memory` argument to `cyto_plot()` to allow interactively positioned label co-ordinates to be passed to `cyto_plot()` when `cyto_plot_save()` in called (as we cannot interactively position labels on these graphics devices).
* Added new `cyto_unmix_compute()` and `cyto_unmix()` to compute and apply spectral unmixing matrices to spectral flow cytometry data.
* `channel_match` argument has now been removed from `cyto_spillover_compute()` and `cyto_spillover_edit()` as CytoExploreR will automatically search in the current working directory for a file called `Compensation-Details` to import these details. 
* Bagwell method in `cyto_spillover_compute()` has been updated to better handle controls from different sources (i.e. beads or cells). `cyto_spillover_compute()` will also automatically handle controls for the same channels split between groups by selecting the control with the highest spillover coefficients. If multiple unstained controls are supplied within a group,  `cyto_spillover_compute()` will use the control with the lowest CVs across all channels. `cyto_spillover_compute()` will now suggest a channel for each control if this information has not been manually supplied. A non-interactive gating method has been added to the Bagwell method to automatically gate negative and positive populations using `openCyto::mindensity()`.
* `cyto_transform()` can now perform inverse transformations on GatingHierarchies and GatingSets by applying transformations to gates as well.
* `cyto_transformers_define()` now accepts custom transformers through a named list of the form `list('PE-A' = list('transform' = asinh, 'inverse' = sinh))` where `asinh` and `sinh` are the names of the tranformer functions.
* `cyto_save()` can now export multiple populations in a single call by writing FCS files to separate folders. `cyto_save()` now calls `cyto_data_extract()` to allow improved flexibility over how the data should be exported (e.g. events, coerce, inverse, copy etc.).
* `cyto_data_extract()` now appropriately handles cytoframes by using the `uri` instead of `cyto_names()`.
* `cyto_load()` and `cyto_setup()` now accept a named list of matrices through the `path` argument to streamline conversion into cytoset objects for use within CytoExploreR. These functions also get new `fixed` and `ignore.case` arguments to provide more flexibility over file selection through the `select` and `exclude` arguments.
* `cyto_plot_contour()` now uses a faster binned kernel density estimator that generates contour lines with ~3.5x improved resolution.
* Added new `cyto_plot_animate()` to create animated GIFs of plots generated by `cyto_plot()`.
* Added new `cyto_ggplot()` to allow conversion of plots created by `cyto_plot()` into grob or ggplot objects, that can easily be arranged with other such objects using `cowplot` or `patchwork` packages.
* Add support for `PeacoQC` anomaly detection algorithm within `cyto_clean()`.
* `cyto_calibrate()` and `cyto_calibrate_reset()` have been renamed to `cyto_plot_calibrate()` and `cyto_plot_calibrate_reset()` to make it easier to find these functions and provide more consistency with other `cyto_plot_()` family members. `cyto_plot_calibrate()` now supports calibration on a per channel basis and users can manually supply channel limits through the new `limits` argument.
* Added new set of `key` arguments to `cyto_plot()` to appropriately annotate point colour scales with counts or fluorescent intensities. This is a major update to `cyto_plot()` which now uses a fixed scale between plots for point colour scales by default (see new `key_scale` argument).
* Added new `cyto_infinity_screen()` function to add support for infinityFlow within CytoExploreR.
* Added new `order` argument to `cyto_plot_profile()` and `cyto_plot_explore()` to allow control over whether plots should be generated by `"channels"` or `"groups"`.
* Major improvements to `cyto_plot_compensation()` to better handle compensated and/or transformed data. Added support for fitting robust linear models and displaying spillover coefficients on plots.
* Added new `"trim"` option to `axes_limits` argument in `cyto_plot()` to trim lowest 1% of events from each channel when computing axes limits.
* Autospill is now the default method used to compute the spillover matrix in `cyto_spillover_compute()`. Users can control which method to use (either Bagwell or Roca) through the new `type` argument.
* Add `copy` argument to `cyto_compensate()` to apply compensation to a copy of the supplied data.
* `cyto_compensate()` now supports spillover matrices in long format (i.e. non-square matrix with 3 columns of form | channel | channel | value |).
* Add new `cyto_plot_new_page()` function to allow plotting on new page when previous page contains empty panels.
* Add optional copy argument to `cyto_transform()` to apply transformations to  copy of the supplied data.
* `cyto_map()` has been re-written to add support for extracting events on a per sample basis. A new `scale` argument has been added to optionally scale each channel to `mean`, `range` or `zscore`. A new `label` argument has been added to allow more fine control over the names of the dimension-reduced parameters.
* BREAKING CHANGE: `display` argument has been replaced with `events` in all functions for consistency.
* `path` argument in `cyto_load()` has been updated to allow support for multiple directories or paths to files.
* A new `overwrite` argument has been added to `cyto_save()` to provide a non-interactive way of controlling how existing directories should be handled when attempting to save files.
* Add new `id` argument to `cyto_split()` to allow more control over how samples are split. The sample IDs can now be stored at the `cytoframe` level or passed manually as a vector or factor to `id`. The names will automatically be inherited from factors where appropriate.
* Add new logical `append` argument to `cyto_markers_extract()` to optionally return extracted markers in the form `<marker> channel`.
* Add new family of `cyto_func_()` functions to more elegantly handle the calling of internal/external functions.
* `cyto_spillover_edit()` now has up/down buttons next to channel selectors to make it easier to scroll through channels.
* `cyto_select()` can now extract samples based on partial string matches for file names. For example, `cyto_select(gs, c("7AAD", "PE"))` will extract all samples with `7AAD` or `PE` in their file names. A new `cyto_match()` function has been added to return the indices of samples that match the selection criteria.
* `cyto_gate_edit()` now appropriately handles the `negate` argument and allows addition/removal of negated gates.
* `group_by` argument in `cyto_gate_draw()` and `cyto_gate_edit()` has been renamed to `merge_by` to maintain consistency with `cyto_plot()`.
* Added new `cyto_gate_extract()` function to facilitate the extraction and formatting of gates extracted from GatingHierarchies, GatingSets or gatingTemplates.
* `cyto_nodes()` family of functions have been updated to include support for extracting and manipulating nodes stored in gatingTemplates. A new `terminal` argument has been added to return the paths of nodes which have no descendants.
* Added `cyto_nodes_ancestor()` to make it easy to get the most recent common ancestor for a collection of populations. 
* `cyto_gate_bool()` now checks to see if the boolean gate is already defined within the gatingTemplate and asks the user if they would like to update the logic for this gate. As an alternative, boolean gates can be modified manually by directly editing the ntries in the gatingTemplate using `cyto_gatingTemplate_edit()`.
* `cyto_gate_draw()` now uses a single method that has support for cytoset and GatingSet objects. The cytoset method simply returns the constructed gates in a list, whilst the GatingSet method adds the gates to the GatingSet, updates the gatingTemplate and returns the updated GatingSet.
* `cyto_clean()` now supports flowAI, flowClean and flowCut.
* Added new `cyto_require()` function to handle imports of external packages not included with CytoExploreR.
* `cyto_plot()` now has a new colourblind friendly palette for 2D scatter plots.
* `cyto_gatingTemplate_create()` now performs a directory check to ensure that the specified gatingTemplate does not already exist.
* Add `active` argument to `cyto_gatingTemplate_create()` to allow users to automatically set a newly created gatingTemplate as the active gatingTemplate.
* `cyto_gatingTemplate_select()` is now defunct in favour of `cyto_gatingTemplate_active()` which can now either set or retrieve the current active gatingTemplate.
* Handling of negated gates within `cyto_gate_draw()` has been significantly improved.
* Substantial speed improvements have been made to `cyto_gate_draw()` as a result of simplified internals and the use of the same data preparation method used within `cyto_plot()`.
* The `group_by` argument in `cyto_gate_draw()` has been renamed to `merge_by` to match `cyto_plot()`. 
* Added new `label_text_col_alpha` argument to `cyto_plot()` to allow control over transparency of label text.
* Added new `cyto_gate_apply()` function to handle all gating operations within `CytoExploreR`. This handy function applies gates to cytoframes or cytosets and returns a list of gated populations.
* cytoframe method for `cyto_save()` has been removed, users should use cytosets instead.
* `cyto_split()` now returns a cytoset and also accepts cytosets containing merged cytoframes.
* `popup` is now set to `TRUE` by default in `cyto_plot()` and support has been added for modifying the size of the pop-up graphics device through the `popup_size` argument. This argument has also been added to `cyto_plot_theme()` so that the default size of the graphics device can be modified globally.
* The entire family of `cyto_transformer()` functions have been completely re-written. The old transformer function have now been merged into a single function called `cyto_transformers_define()`. This function replaces the old `cyto_transformer_log()`, `cyto_transformer_arcsinh()`, `cyto_transformer_biex()` and `cyto_transformer_logicle()` functions which are now defunct. The old `cyto_transformer_combine()` and `cyto_transformer_extract()` functions have also been renamed to instead use the more consistent `cyto_transformers` prefix.
* `cyto_gate_remove()` now handles removal of populations defined by quadrant gates.
* All density arguments within `cyto_plot()` have been renamed to instead carry a `hist` prefix. A new `hist_stat` argument has been added to control the statistic to display on the plot, options include `count`, `density` or `percent`. Data outside plot limits is now excluded from the kernel density computation to improve visualisation of the data.
* Significant improvements have been made to the handling of 2D density colours within `cyto_plot()`. The resolution has been improved by increasing the number of bins used to construct the grid when computing the binned 2D density kernel (new default is set to 200 bins). The grid is also now increased proportionately if there is data outside the plot limits, thus prevent stretching of bin colours within the plot (preventing a box-like appearance).
* `cyto_plot_contour()` has been re-written to accept any cytometry data type and to restrict 2D density computation to be within plot limits. 
* `cyto_plot()` gains significant speed improvements due to changes in the way that data is sampled and merged internally. A new function `cyto_sample_n()` is now used compute sample sizes for cytoframes prior to merging.
* Fix labelling of transformed axes for values less than zero.
* `cyto_sample()` has been updated to add support for individually sampling cytoframes within a cytoset.
* `cyto_merge_by` has gained a new `format` argument to control how grouped data is returned. 
* `cyto_load()` now barcodes each event with an ID so that `cyto_plot()` can plot the exact same events when using `overlay`. This behavior is different to `flowWorkspace::load_cytoset_from_fcs()` which means that users will need to use the CytoExploreR `cyto_load()` or `cyto_setup()` APIs instead.
* `cyto_barcode()` has been re-written to allow replacement of existing barcodes.
* Directories are now excluded from the files read in by `cyto_load()`.
* `cyto_stats_compute()` also has enhanced statistical capabilities with support added for standard deviation (SD), robust standard deviation (rSD), coefficient of variation (CV), robust coefficient of variation (rCV), quantile, area under curve (AUC) and range.
* `cyto_stats_compute()` now accepts custom functions and experimental details can be optionally excluded from the output. Statistics that rely on density distributions now use the same bandwidth for each cytoframe.
* Rename all `density` arguments to instead use `hist` as a prefix.
* Add legend fro point colour scale to `cyto_plot()`.
* Add support for adding grid lines to `cyto_plot()` through `grid` argument and customisation through `grid_line_width`, `grid_line_col` and `grid_line_alpha` arguments.
* Add new `cyto_plot_par()` function to handle graphical parameters within `cyto_plot()`. Refine behaviour of all `cyto_plot()` helper functions to ensure graphics devices are opened when expected. `cyto_plot_layout()` is now defunct, `layout` can now be directly passed to `cyto_plot_new()`, `cyto_plot_custom()` or `cyto_plot_complete()`.
* Support for scattermore is now integrated in `cyto_plot_point()` through the `point_fast` argument, which is set to FALSE by default.
* Add new `cyto_gate_copy()` function to support copying of gates from one node to another.
* Add support for partial matches in `cyto_channels_extract()` and `cyto_markers_extract()`. This means that partial matching of channels/markers is now supported by `cyto_plot()` as well.
* Add support for extracting channels from a list of cytometry objects in `cyto_channels()` and `cyto_markers()`. Add new replacement method for `cyto_channels()`.
* Remove old data editor in favour of `DataEditR::data_edit()`.
* Add new `cyto_cbind()` function to handle addition of new parameters to cytometry objects (retaining support for flowFrames/flowSets).
* `cyto_apply()` has been re-written to accept different data inputs (e.g. matrix, cytoframe, cytoset, column or row) and the output is more intuitively formatted. A new `slot` argument has been added to allow extraction of data from a specific slot of a function output.
* `cyto_calibrate()` now uses quantile calibration by default and the upper quantile limits has been reduced from 0.99 to 0.95.
* Barcoding of GatingSets is now supported by `cyto_barcode()`.
* `cyto_beads_sample()` is now defunct in favour of the more flexible `cyto_sample_to_node()` function.
* Add new `cyto_sort_by()` function to sort samples based on experimental variables.
* Add new `cyto_groups()` function to return the names or details of experimental groups.
* `cyto_extract()` is now defunct in favour of the new `cyto_data_extract()` which has for extracting multiple populations. New arguments `format`, `split`, `coerce`, `markers`, `trans`, `inverse` and `sample` provide flexibility in how the extracted population should be formatted.
* `cyto_names()` has been updated to extract sample information at the level of cytosets only, attempts to extract this information from flowFrames/cytoframes will generate an error.
* `cyto_details()` now has additional `convert`, `drop` and `factor` arguments to allow easy formatting of the returned data.frame. Support for extracting `cyto_details()` from flowFrame/cytoframe objects has been removed.
* `cyto_transform()` now appropriately handles transforming subsetted data where some parameters defined in the transformerList may be missing in the data.
* `cyto_details()` called on a cytoframe now returns a data.frame to match the behaviour of cytosets, GatingHierarchies and GatingSets.
* Add new `cyto_class()` to get or check the class of cytometry objects.
* Remove directories from file paths in `cyto_load()`. Re-factor code to ensure loading of saved GatingSets.

# CytoExploreR 1.1.0

* Fix `cyto_gate_bool()` to anchor new nodes to most recent ancestor in the GatingSet.
* Fix labelling of transformed axes within `cyto_plot()`.
* Fix statistics computation in `cyto_plot()` when gates and overlays are present.
* Fix `cyto_plot_gating_scheme()` to allow the use of custom colour scale when `back_gate = TRUE`.
* Make sure `cyto_load()` is able to read in saved GatingSets.

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