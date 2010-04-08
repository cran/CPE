.First.lib <- function(lib, pkg) {
   library.dynam("CPE", pkg, lib)
   cat("CPE 1.4.1 loaded\n")
}

