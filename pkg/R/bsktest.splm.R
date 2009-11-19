`bsktest.splm` <-
function(x, w, index=NULL, test=c("CLMlambda","CLMmu"), ...){
	
	switch(match.arg(test), CLMlambda = {

    bsk = clmltest.model(x,w, index)

  }, CLMmu = {

    bsk = clmmtest.model(x,w, index )

  })

  return(bsk)
}

