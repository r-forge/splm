`conclikpan` <-
function(lambda, eig, e0e0, e1e1,e0e1, N,T,NT,quiet=quiet)
{
	Nsig <- e0e0 - 2*lambda*e0e1 + lambda*lambda*e1e1
	sigma2 <- Nsig/NT
	if (is.complex(eig)) det <-  Re(prod(1 - lambda*eig)) 
	else det <-  prod(1 - lambda*eig)
	ret <-  T * log(det) - ( (NT)/2 )*log(2*pi*sigma2) - (1/(2*sigma2))*Nsig
	#ret <- - (NT/2)*log(Nsig)  + T * log(det)  
#print(T*log(det))
	  if (!quiet) 
        cat("(eigen) rho:\t", lambda, "\tfunction value:\t", ret, 
            "\n")

	ret
}

