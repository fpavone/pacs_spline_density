.onAttach <- function(libname, pkgname)
  packageStartupMessage("splineDensity 0.0.0.23 loaded\nCopyright Noi 2018")

.onUnload <- function(libpath)
  library.dynam.unload("splineDensity",  libpath)