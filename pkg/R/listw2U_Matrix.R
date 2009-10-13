`listw2U_Matrix` <-
function(lw){
as(as(0.5 * (lw + t(lw)), "symmetricMatrix"), "CsparseMatrix")
}

