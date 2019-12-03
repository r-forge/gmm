#####  All S4 classes of the package are defined here
######################################################


## Causal Model Classes

setClass("causalGel", contains="functionGel")

setClass("causalData", representation(momType="character",
                                      balCov="character",
                                      balMom="numericORNULL",
                                      ACTmom="integer",
                                      reg="data.frame",
                                      bal="data.frame"))

setClass("causalGelfit", contains="gelfit")


## converters

