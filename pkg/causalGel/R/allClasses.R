#####  All S4 classes of the package are defined here
######################################################


## Causal Model Classes

setClass("causalModel", contains="functionModel")

setClass("rcausalModel", contains="rfunctionModel")

setClass("causalData", representation(momType="character",
                                      balCov="character",
                                      balMom="numericORNULL",
                                      ACTmom="integer",
                                      reg="data.frame",
                                      bal="data.frame"))

setClass("causalGelfit", contains="gelfit")

setClass("causalfit", representation(estim="numeric",
                                     se = "numeric",
                                     coefNames = "character",
                                     type="character",
                                     method="character",
                                     form="list",
                                     details="list",
                                     info="list",
                                     data="data.frame",
                                     call="callORNULL"))

setClass("causalGelfit", contains="gelfit")

setClass("summaryCausalfit", representation(coef="matrix",
                                            type="character",
                                            method="character",
                                            form="list",
                                            details="list",
                                            info="list"))

## converters

setAs("rcausalModel", "causalModel",
      function(from) {
          obj <- as(from, "momentModel")
          new("causalModel", obj)})


