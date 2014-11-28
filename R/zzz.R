
.onLoad <- function(libname, pkgname) {
    library.dynam("bigtabulate", pkgname, libname)
}

#.noGenerics <- TRUE

.onUnload <- function(libpath) {
    library.dynam.unload("bigtabulate", libpath)
}
