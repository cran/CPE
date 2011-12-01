.First.lib <- function(lib, pkg) {
   library.dynam("CPE", pkg, lib)
   cat("CPE 1.4.2 loaded\n")
}

