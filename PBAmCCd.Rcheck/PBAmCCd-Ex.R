pkgname <- "PBAmCCd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('PBAmCCd')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("g")
### * g

flush(stderr()); flush(stdout())

### Name: g
### Title: Log of the Inverse Additive Log-ratio
### Aliases: g

### ** Examples

vec <- sample(1:250, 9)
ref <- vec[length(vec)]
alr.vec <- log(vec[-length(vec)]/ref)

g(alr.vec)




cleanEx()
nameEx("rmse")
### * rmse

flush(stderr()); flush(stdout())

### Name: rmse
### Title: Root Mean Square Error
### Aliases: rmse

### ** Examples

x <- sample(1:20, size = 4)
y <- sample(1:20, size = 4)

rmse(x,y)




cleanEx()
nameEx("rmse_by_row")
### * rmse_by_row

flush(stderr()); flush(stdout())

### Name: rmse_by_row
### Title: Row Root Mean Square Error
### Aliases: rmse_by_row

### ** Examples

x <- matrix(sample(1:500, size = 21), nrow=3)
y <- matrix(sample(1:500, size = 21), nrow=3)

rmse(x,y)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
