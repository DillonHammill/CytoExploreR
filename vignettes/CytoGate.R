## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", message = FALSE, warning = FALSE)

## ---- echo = TRUE--------------------------------------------------------
library(openCyto)
library(CytoGate)

registerPlugins(fun = gate_draw, methodName = "DrawGate")
listgtMethods()


## ---- echo = TRUE--------------------------------------------------------
fs <- Activation

gs <- GatingSet(fs)

## ---- echo = TRUE--------------------------------------------------------
spill <- fs[[1]]@description$SPILL
gs <- compensate(gs, spill)

channels <- colnames(fs)
trans.channels <- channels[!channels %in% c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","Time")]

trans <- estimateLogicle(gs[[1]], trans.channels)
gs <- transform(gs, trans)

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "Cells", parent = "root", pop = "+", dims = "FSC-A,SSC-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='polygon'", collapseDataForGating = TRUE, groupBy = 8
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "PE+", parent = "root", pop = "+", dims = "FSC-A,SSC-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='rectangle'", collapseDataForGating = TRUE, groupBy = 2
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "PE+", parent = "root", pop = "+", dims = "PE-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='interval'", collapseDataForGating = TRUE, groupBy = 2
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "PE+", parent = "root", pop = "+", dims = "APC-Cy7-A,PE-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='interval',axis='y'", collapseDataForGating = TRUE, groupBy = 2
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "PE+", parent = "root", pop = "+", dims = "PE-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='threshold'", collapseDataForGating = TRUE, groupBy = 2
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "PE+FITC+", parent = "root", pop = "+", dims = "Alexa Fluor 488-A,PE-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='threshold'", collapseDataForGating = TRUE, groupBy = 2
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "PE+FITC+", parent = "root", pop = "+", dims = "Alexa Fluor 488-A,PE-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='ellipse'", collapseDataForGating = TRUE, groupBy = 2
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "DN,FITC+,AF700+FITC+,AF700+", parent = "root", pop = "*", dims = "Alexa Fluor 488-A,Alexa Fluor 700-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='quadrant'", collapseDataForGating = TRUE, groupBy = 2
#  )

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  template <- add_pop(
#    gs, alias = "PE+,FITC+", parent = "root", pop = "*", dims = "Alexa Fluor 488-A,PE-A", gating_method = "DrawGate",
#    gating_args = "subSample=25000,gate_type='polygon',N=2", collapseDataForGating = TRUE, groupBy = 2
#  )

