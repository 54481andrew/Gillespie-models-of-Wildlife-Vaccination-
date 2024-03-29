
#h <- function(t,parms){with(parms, b0/(1 + ((t%%T)/tb)^40))}
h <- function(t,parms){with(parms, b0*(t%%T < tb) )}

rhs.fun = function(t,y,parms){
    S <- y[1]
    Iv <- y[2]
    Ip <- y[3]
    V <- y[4]
    P <- y[5]
    with(parms,{
        dS = h(t,parms) - Bp*S*Ip - d*S
        dIv = -(gamv + d)*Iv - Bp*Iv*Ip
        dIp = Bp*S*Ip + Bp*Iv*Ip - (gamp + d)*Ip
        dV = gamv*Iv - d*V
        dP = gamp*Ip - d*P
        
        return(list(c(dS,dIv,dIp,dV,dP)))
    })
}

rhs.freq.fun = function(t,y,parms){
    S <- y[1]
    Iv <- y[2]
    Ip <- y[3]
    V <- y[4]
    P <- y[5]
    N <- sum(y)
    with(parms,{
        dS = h(t,parms) - Bp*S*Ip/N - d*S
        dIv = -(gamv + d)*Iv - Bp*Iv*Ip/N
        dIp = (Bp*S*Ip + Bp*Iv*Ip)/N - (gamp + d)*Ip
        dV = gamv*Iv - d*V
        dP = gamp*Ip - d*P
        
        return(list(c(dS,dIv,dIp,dV,dP)))
    })
}


vaccinate <- function(t,y,parms){
    S <- y[1]
    Iv <- y[2]
    Ip <- y[3]
    V <- y[4]
    P <- y[5]
    

    with(parms,{
	nvacc <- round(min(S, Nv*S/(S + Iv + Ip + V + P)))
	S   <- S - nvacc
	Iv  <- Iv + nvacc
	y <- c(S, Iv, Ip, V, P)
	return(y)
	})
}



