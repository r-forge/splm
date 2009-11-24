`bsktest.lm` <-
function(x, w, index=NULL, test=c("SLM1","SLM2","LMJOINT"), ...){
	
	switch(match.arg(test), SLM1 = {

    bsk = slm1test.model(x,w, index)

  }, SLM2 = {

    bsk = slm2test.model(x,w, index )

  }, LMJOINT = {

    bsk = LMHtest.model(x,w, index)

  })

  return(bsk)
}

