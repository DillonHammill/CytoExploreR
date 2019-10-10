# CytoRSuite 1.0.0

* Improved handling of lower axes limits to display negative values.
* Improved support for handling samples with few events.
* Updated `gate_draw`, `gate_edit` and `gate_remove` to allow easier switching between these functions. Users can now replace `draw` in their `gate_draw` code with `edit` or `remove` to modify drawn gates without any additional code.
* Supplying an empty character string to the `alias` argument of `cyto_plot` (e.g. alias = "") will automatically plot all gates constructed in the supplied channels. This way users don't have to specify each gate by name.
* Added support for spillover matrices has been added through `spillover_spread_compute` which utilises a similar API to `spillover_compute`.
* Modified colour scheme in `spillover_edit` to improve visibility.
* Improved layout for `cyto_plot_gating_scheme` when there isn't a `legend`.
* Updated `cyto_plot` to allow plotting of all 2-D gate objects in a single dimension. The minimum and maximum gate co-ordinates in the supplied channel will be used to construct a 1-D rectangleGate for plotting.
* Labels are now adjusted to prevent overlap in `cyto_plot_label` and `cyto_plot` when multiple gates are supplied.
* `gate_draw` now restricts large flowSet or GatingSet objects to 20 random samples to improve processing speed.
* New function `cyto_plot_save` provides a simpler method for saving high resolution images. `cyto_plot_save` should be called prior to plotting.
* `gate_rename` provides an easy way to update gate names in the GatingHierarchy/GatingSet and associated gatingTemplate.
* `cyto_stats_compute` has been completely revamped to improve speed and to return tidyverse-friendly tibbles in either wide or long formats.

* cyto_plot 1D - center legend, return smoothing functionality, no longer require modal normalisation for overlays, now adds empty space for samples with missing data (NaN%), fix behavior of titles when overlays are present.
* interactively position labels :)

# CytoRSuite 0.9.9

* `cyto_plot_label` now accepts stat `freq` instead of `percent` for consistency with `cyto_stats_compute`.
* `density_layers` argument for `cyto_plot` is now operational. This argument designates the number of layers to include per plot for stacked density distributions. 
* `cyto_markers` has been updated to allow editing of channel names as well.
* `cyto_stats_compute` now calculates the robust CV for stat = "CV".
* New function `gatingTemplate_edit` to facilitate interactive editing of the gatingTemplate. The feature is designed for adding boolean, reference or automated gates to existing gatingTemplates.
* Improved support for editing gates constructed with multiple grouping variables using `gate_edit`.
* `cyto_annotate` to save pData information to a csv file ("Experiment-Details.csv") for future use.
* `cyto_annotate` can now accept a csv file containing experiment details instead of manually re-entering this information.
* `cyto_markers` updated to save markers to csv file ("Experiment-Markers.csv") for future use.
* Added Contributors Code of Conduct.
* New 'Manual & Automated Gating' vignette.

# CytoRSuite 0.9.8

* Completely updated CytoRSuite to comply with rOpenSci naming requirements.
* Users can now draw gates per sample or group using `group_by` in `gate_draw`. To facilitate this change
major changes were made to the way the gatingTemplate is handled and manipulated in the package. Users of 0.9.5 should use `gatingTemplate_convert` to convert their old gatingTemplates to this new format.
* Full support for back-gating is now supported in `cyto_plot` and `cyto_plot_gating_scheme` through `overlay` and `back_gate` arguments.
* New gate tracking feature added to `cyto_plot_gating_scheme` to colour gates and plot borders based on population.
* `spillover_edit` has been updated to include density histograms with statistics for unstained and stained controls. Statistics automatically update as spillover matrix is edited.
* `cyto_plot_compensation` is now incorporated directly into `spillover_edit` in the plots tab.
* `cyto_markers` provides an interactive editable table to assign marker names to channels.
* `cyto_annotate` provides an interactive editable table to annotate samples with experimental details.
* New `cyto_plot` functions have been added to layer overlays (`cyto_plot_overlay`), gates (`cyto_plot_gate`), labels (`cyto_plot_label`) or contour lines (`cyto_plot_contour`) onto existing plots.
* Many more changes to improve the stability of the package. Refer to the documentation for additional changes.

# CytoRSuite 0.9.5

* Added a `NEWS.md` file to track changes to the package.