.First.lib <- function(libname, pkgname, where) {


  where <- match(paste("package:", pkgname, sep=""), search())

  require(Biobase, quietly=TRUE) ##Biobase uses methods

    if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addPDF2Vig("affycomp")
    }



  cacheMetaData(as.environment(where))

}
