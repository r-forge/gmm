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


## converters

setAs("rcausalModel", "causalModel",
      function(from) {
          obj <- as(from, "momentModel")
          new("causalModel", obj)})


