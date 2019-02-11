# CytoRSuite 0.9.9

* Improved support for editing constructed with multiple grouping variables using `gate_edit`.
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
