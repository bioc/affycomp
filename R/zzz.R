.onLoad <- function(libname, pkgname) require("methods")

.onAttach <- function(libname, pkgname) {

    if (.Platform$OS.type == "windows" && interactive()
        && .Platform$GUI == "Rgui") {
        addVigs2WinMenu("affycomp")
    }
}
