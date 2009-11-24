`bsktest.formula` <-
function(x, data, w, test=c("SLM1","SLM2","LMJOINT","CLMlambda","CLMmu"), index=NULL, ...){
  

switch(match.arg(test), SLM1 = {

    bsk = slm1test(x, data, index,  w)

  }, SLM2 = {

    bsk = slm2test(x, data, index,  w)

  }, LMJOINT = {

    bsk = LMHtest(x, data, index,  w)

  }, CLMlambda = {

    bsk = clmltest(x, data, index,  w)

  }, CLMmu = {

    bsk = clmmtest(x, data, index,  w)

  })

  return(bsk)

}

