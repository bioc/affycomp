.First.lib <- function(libname, pkgname, where) {
  

  where <- match(paste("package:", pkgname, sep=""), search())

  require(Biobase, quietly=TRUE) ##Biobase uses methods

  cacheMetaData(as.environment(where))

}
