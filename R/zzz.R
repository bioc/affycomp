.onAttach <- function(libname, pkgname) {
    if (.Platform$OS.type == "windows" && interactive()
        && .Platform$GUI == "Rgui") {
        addVigs2WinMenu("affycomp")
    }
}

## Fixing 'no visible' notes
.myDataEnv <- new.env(parent=emptyenv())

getData <- function(dataset, package=NULL){
  if (!exists(dataset, .myDataEnv))
    data(list=dataset, package=package, envir=.myDataEnv)
  .myDataEnv[[dataset]]
}
