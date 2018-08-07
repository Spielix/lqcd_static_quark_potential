# Header for command line arguments:
# #!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#######################################

lm.avail <- require(minpack.lm)
parallel.avail <- require(parallel)

if(!lm.avail){
	stop("minpack.lm not found!\n", call.=FALSE)
}

## chi function
mychi <- function(f, par, x, y, dy) {
  (y - f(par=par, x=x))/dy
}

## fits, returns fit-results and chi_sq
wrapper <- function(y, ...) {
  res <- nls.lm(fn=mychi, y=y, ...)
  return (c(res$par, res$rsstrace[length(res$rsstrace)]))
}

## define fit-routine
jfit <- function(x, y, dy, f, par, boot.R=1500, show=TRUE){
	bootstrapsamples <- array(NA, dim=c(boot.R+1, length(y)))
	bootstrapsamples[1,] <- y
	rr <- c(2:(boot.R+1))
	if(parallel.avail){
		bootstrapsamples[rr,] <- mcmapply(rnorm, n=boot.R, mean = y, sd = dy)
	}else{
		bootstrapsamples[rr,] <- mapply(rnorm, n=boot.R, mean = y, sd = dy)
	}
	boot.res <- apply(bootstrapsamples, MARGIN=1, FUN=wrapper, x=x, dy=dy, par=par, f=f)
	par.booterrors <- apply(boot.res, 1, sd)
	par.bootmean <- apply(boot.res, 1, mean)
	if(show){
		print("Parameters..., chi_sq")
		print(boot.res[,1])
		print(par.booterrors[c(1:(length(par.booterrors)-1))])
		par.bootcor <- cor(t(boot.res[c(1:(length(par))),]), t(boot.res[c(1:(length(par))),]))
		print("")
		print("Corralation matrix")
		print(par.bootcor)
	}
	## res[1,] = fit results, chi_sq
	## res[2,] = fit errors, std dev of chi_sq
	## res[3,] = mean of the bootstrap results and chi_sq
	res <- array(NA, dim=c(4, length(par)+1))
	res[1,] <- boot.res[,1]
	res[2,] <- par.booterrors
	res[3,] <- par.bootmean
	res[4,] <- par.bootcor[2,3]
	return (res)
}

jplot <- function(x, y, dy, dx, col="black", ...) {
	plot(x=x, y=y, col=col, ...)
	arrows(x0=x, y0=y-dy, x1=x, y1=y+dy, length=0.01, angle=90, code=3, col=col)
	if(!missing(dx)){
		arrows(x0=x-dx, y0=y, x1=x+dx, y1=y, length=0.01, angle=90, code=3, col=col)
	}
}

jlines <- function(f, par, samples=400, data, xmin, xmax, ...){
	if(missing(data)){
		x = seq(xmin, xmax, length.out=samples)
	}else{
		x = seq(min(data), max(data), length.out=samples)
	}
	y = f(par, x)
	lines(x, y, ...)
}

##############################################################################
# Global definitions end
# Local script starts
##############################################################################

# Example for usage
energy <- function(par, x){
	par[1]+par[2]*x+par[3]/x
}

data = read.table(args[1])

par = rep(1,3)
res <- jfit(x=data[[1]], y=data[[2]], dy=data[[3]], f=energy, par=par)
a = sqrt(res[1,2] / (1.65 + res[1,3]))/2.
aerr = sqrt((res[2,2]/(8.*a*(1.65 + res[1,3])))^2+(res[2,3]*a/(2*(1.65+res[1,3])))^2+res[4,1]*res[2,2]*res[2,3]*(-1./(16.*(1.65+res[1,3])^2)))
temperatur = 1.97e-13/(4.*a*1e-15)
temperaturerr = temperatur*aerr/a
print(a)
print(aerr)
print(temperatur)
print(temperaturerr)
jplot(data[[1]], data[[2]], data[[3]], xlab="m_s", ylab="E")
jlines(energy, res[1,c(1:length(par))], data=data[[1]])
