#####  All S4 classes of the package are defined here
######################################################


## Causal Model Classes

setClass("causalGel", contains="functionGel")

setClass("rcausalGel", contains="rfunctionGel")

setClass("causalData", representation(momType="character",
                                      balCov="character",
                                      balMom="numericORNULL",
                                      ACTmom="integer",
                                      reg="data.frame",
                                      bal="data.frame"))

setClass("causalGelfit", contains="gelfit")


## converters

setAs("rcausalGel", "rgmmModels",
      function(from) {
          as(as(from, "rgelModels"), "rgmmModels")})

setAs("rcausalGel", "causalGel",
      function(from) {
          obj <- as(from, "gelModels")
          new("causalGel", obj)})


