.onLoad <- function(libname, pkgname) require("methods")

.onAttach <- function(libname, pkgname) {

    if (.Platform$OS.type == "windows" && interactive()
        && .Platform$GUI == "Rgui") {
        addVigs2WinMenu("affycomp")
    }
}

## Fixing 'no visible' notes
.myDataEnv <- new.env(parent=emptyenv())

isLoaded <- function(dataset)
  exists(dataset, .myDataEnv)

getData <- function(dataset){
  if (!isLoaded(dataset))
    data(list=dataset, envir=.myDataEnv)
  .myDataEnv[[dataset]]
}
