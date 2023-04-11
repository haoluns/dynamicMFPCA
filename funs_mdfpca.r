

library(dplyr)
library(fda)


fd2d = function(betas,basis2d){
	
	prod1 = inprod(basis2d$basis1, basis2d$basis1,0,0)
	prod2 = inprod(basis2d$basis2, basis2d$basis2,0,0)
	prod0 = inprod(basis2d$basis0, basis2d$basis0,0,0)
	
	return(list(betas=betas,basis1=basis2d$basis1,basis2=basis2d$basis2,basis0=basis2d$basis0,prod1=prod1,prod2=prod2,prod0=prod0))
}


fd2d = function(betas,basis1,basis2,basis0){
	
	prod1 = inprod(basis1, basis1,0,0)
	prod2 = inprod(basis2, basis2,0,0)
	prod0 = inprod(basis0, basis0,0,0)
	
	return(list(betas=betas,basis1=basis1,basis2=basis2,basis0=basis0,prod0=prod0,prod1=prod1,prod2=prod2))
}


create.2d.basis = function(b1,b2){

	return (list(basis1=b1,basis2=b2))
}

norm.fd2d=function(fd2dobject){

betas = fd2dobject$betas
basis1 = fd2dobject$basis1
basis2 = fd2dobject$basis2
basis0 = fd2dobject$basis0

prod0 = fd2dobject$prod0
prod1 = fd2dobject$prod1
prod2 = fd2dobject$prod2

numbasis0=basis0$nbasis
numbasis1=basis1$nbasis
numbasis2=basis2$nbasis



temp=0
for(k in 1:numbasis0){
 for(i in 1:numbasis1){
 for(j in 1:numbasis2){
 for(kp in 1:numbasis0){
 for(ip in 1:numbasis1){
 for(jp in 1:numbasis2){
 
 temp = temp +betas[i,j,k]*betas[ip,jp,kp]*prod1[i,ip]*prod2[j,jp]*prod0[k,kp]
}
 }
 }
 }
 }
 }
 
return (temp)
}





inprod.fd2d=function(fd2dobject1,fd2dobject2){

betas1 = fd2dobject1$betas

basis0 = fd2dobject1$basis0
basis1 = fd2dobject1$basis1
basis2 = fd2dobject1$basis2

prod0 = fd2dobject1$prod0
prod1 = fd2dobject1$prod1
prod2 = fd2dobject1$prod2

betas2 = fd2dobject2$betas

numbasis0=basis0$nbasis
numbasis1=basis1$nbasis
numbasis2=basis2$nbasis



temp=0
 for(k in 1:numbasis0){
 for(i in 1:numbasis1){
 for(j in 1:numbasis2){
 for(kp in 1:numbasis0){
 for(ip in 1:numbasis1){
 for(jp in 1:numbasis2){
  
  temp = temp +betas1[i,j,k]*betas2[ip,jp,kp]*prod1[i,ip]*prod2[j,jp]*prod0[k,kp]
}
 }
 }
 }
 }
 }
return (temp)
}



inprod.fd2dx=function(fd2dobject1,fd2dobject2){

betas1 = fd2dobject1$betas

basis0 = fd2dobject1$basis0
basis1 = fd2dobject1$basis1
basis2 = fd2dobject1$basis2

prod0 = fd2dobject1$prod0


basis01 = fd2dobject1$basis0
basis11 = fd2dobject1$basis1
basis21 = fd2dobject1$basis2

basis02 = fd2dobject2$basis0
basis12 = fd2dobject2$basis1
basis22 = fd2dobject2$basis2

prod0 = inprod(basis01, basis02,0,0)
prod1 = inprod(basis11, basis12,0,0)
prod2 = inprod(basis21, basis22,0,0)


betas2 = fd2dobject2$betas

numbasis0=basis0$nbasis
numbasis1=basis1$nbasis
numbasis2=basis2$nbasis



temp=0
 for(k in 1:numbasis0){
 for(i in 1:numbasis1){
 for(j in 1:numbasis2){
 for(kp in 1:numbasis0){
 for(ip in 1:numbasis1){
 for(jp in 1:numbasis2){
  
  temp = temp +betas1[i,j,k]*betas2[ip,jp,kp]*prod1[i,ip]*prod2[j,jp]*prod0[k,kp]
}
 }
 }
 }
 }
 }
return (temp)
}





constraint.fd2d=function(fd2dobject1){

betas1 = fd2dobject1$betas

basis0 = fd2dobject1$basis0
basis1 = fd2dobject1$basis1
basis2 = fd2dobject1$basis2


prod0 = fd2dobject1$prod0
prod1 = fd2dobject1$prod1
prod2 = fd2dobject1$prod2

betas2 = fd2dobject1$betas

numbasis0=basis0$nbasis
numbasis1=basis1$nbasis
numbasis2=basis2$nbasis


ct = array(0,c(numbasis1,numbasis2,numbasis0))

temp=0
 for(k in 1:numbasis0){
 for(i in 1:numbasis1){
 for(j in 1:numbasis2){
 
 
 for(kp in 1:numbasis0){
 for(ip in 1:numbasis1){
 for(jp in 1:numbasis2){
  
  ct[i,j,k] = ct[i,j,k] +betas2[ip,jp,kp]*prod1[i,ip]*prod2[j,jp]*prod0[k,kp]
}
 }
 }
 
 }
 }
 }
 
 
return (betatovec(ct))
}



eval.fd2d=function(timei1,timei2,timei0, fd2dobject){


betas = fd2dobject$betas
basis1 = fd2dobject$basis1
basis2 = fd2dobject$basis2
basis0 = fd2dobject$basis0

numbasis0=basis0$nbasis
numbasis1=basis1$nbasis
numbasis2=basis2$nbasis

T = length(timei1)
value=sapply(1:T, function(t){
bv1 =  eval.basis(timei1[t],basis1)
bv2 = eval.basis(timei2[t],basis2)
bv0 = eval.basis(timei0[t],basis0)
temp = 0

for (k in 1:numbasis0){
for (i in 1:numbasis1){
for (j in 1:numbasis2){
	temp = temp +(bv0[k]) * (bv1[i]) * (bv2[j]) * betas[i,j,k]
}}}
return(temp)
})


return (value)

}

RESBREAK=TRUE




betatovec <- function(betas){
nbasis1 = dim(betas)[1]
nbasis2 = dim(betas)[2]
nbasis0 = dim(betas)[3]

beta1 = array(1,c(nbasis1,nbasis2,nbasis0))
o = seq(1,nbasis1*nbasis2*nbasis0)
for(ind in seq(1,nbasis1*nbasis2*nbasis0)){

index = betatovec0(betas,ind,nbasis1,nbasis2,nbasis0)
o[ind] = betas[index$i,index$j,index$k]
}
return(o)
}		



betatovec0<-function(betas,ind,nbasis1,nbasis2,nbasis0){
	
	
	lenv = length(c(betas))
	

	
	K = nbasis0	
	
	
	
	k = ceiling(ind / (lenv/K))
	
	I = nbasis1
	J = nbasis2
	
	i =ceiling((ind - (k-1)*(lenv/K))/J)
	
	j = (ind - (k-1)*(lenv/K)) - (i-1)*J
	
	return(list(i=i,j=j,k=k))
}


eval.image.process<-function(pc_fit,timepoints0,timepoints1,timepoints2){
	o = sapply(1:length(timepoints0),function(i){
		return( eval.fd2d(timepoints1[i],timepoints2[i],timepoints0[i], pc_fit))
	})
	return(o)
}

vectobeta <- function(beta1vec){
beta1 = array(1,c(nbasis1,nbasis2,nbasis0))
for(i in 1:nbasis1){
for(j in 1:nbasis2){
for (k in 1:nbasis0){
	beta1[i,j,k] = vectobeta0(beta1vec,i,j,k)
}
}
}
return(beta1)
}		

vectobeta0<-function(betalongvec,i,j,k){
	
	K = nbasis0
	
	lenv = length(betalongvec)		
	
	vec1 = betalongvec[seq(1+(lenv/K)*(k-1),1+(lenv/K)*(k-1)+(lenv/K)-1)]
	
	
	I = nbasis1
	
	lenv = length(vec1) 
	vec2 = vec1[seq(1+(lenv/I)*(i-1),1+(lenv/I)*(i-1)+(lenv/I)-1)]
	
	
	J = nbasis2
	
	return(vec2[j])
}



library(parallel)
mlapply<-function(a,f){
x =mclapply(a,f,mc.cores=60)
return(x)
}



findmean = function(observed,timepoints1,timepoints2,timepoints0,basis1,basis2,basis0,threshold=1e-6,minit=3){

yem = do.call(c,observed)

xalpha = mlapply(1:length(observed),function(subj_index){
	
	timei1 =timepoints1[[subj_index]]
	timei2 =timepoints2[[subj_index]]
	timei0 =timepoints0[[subj_index]]
	
	xmati  <- lapply(1:length(timei1),function(t){	
	evalb1 = eval.basis(timei1[t], basis1)
	evalb2 = eval.basis(timei2[t], basis2)
	evalb0 = eval.basis(timei0[t], basis0)
	
	
	v = c(evalb0 %x% evalb1 %x% evalb2 )
	return(v)
	})%>%do.call(rbind,.)			
	return(xmati)
		
})%>%do.call(rbind,.)


yem = as.matrix(yem)
library(bigmemory)
library(bigalgebra)
A = xalpha
tA =t(A)
tyem = t(yem)
tyem = as.big.matrix(tyem,type="double")
yem = as.big.matrix(yem,type="double")
A <- as.big.matrix(A,type="double")
tA <- as.big.matrix(tA,type="double")	

M1 = as.matrix(tA %*% A )
M2 = as.matrix(tA %*% yem)

beta1vec = ginv(M1) %*% M2


nbasis0=basis0$nbasis
nbasis1=basis1$nbasis
nbasis2=basis2$nbasis

beta1 = array(0,c(basis1$nbasis,basis2$nbasis,basis0$nbasis))
for(i in 1:nbasis1){
for(j in 1:nbasis2){
for (k in 1:nbasis0){
	beta1[i,j,k] = vectobeta0(beta1vec,i,j,k)
}
}
}


pc_fit = fd2d(beta1, basis1,basis2,basis0)


return(list(beta=beta1,pc_fit = pc_fit))

}


first_FPC_2d_image = function(beta1,observed,timepoints1,timepoints2,timepoints0,basis1,basis2,basis0,threshold=1e-6,minit=3){
	thresh = 1
	it = 1

	value = -1e14
	pc_fit = fd2d(beta1, basis1,basis2,basis0)
	pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
	
	beta1 = pc_fit$betas
	
	
		
	library(parallel)
	mlapply<-function(a,f){
	x =mclapply(a,f,mc.cores=60)
	return(x)
	}

	
		
	while (thresh >  threshold | it<minit){
		beta1vec_before = betatovec(beta1)
		beta1_before = beta1
		value_before = value
		
		iterationxmatlist=mlapply(1:length(observed),function(subj_index){

		timei1 =timepoints1[[subj_index]]
		timei2 =timepoints2[[subj_index]]
		timei0 =timepoints0[[subj_index]]
			
		iterationxmat = eval.fd2d(timei1,timei2,timei0,pc_fit)
		return(iterationxmat)
		})

	
		alpha_fit = function(subj_index){

			
			xmati = iterationxmatlist[[subj_index]]
			model = lm.fit(x=as.matrix(xmati),y=as.matrix(observed[[subj_index]]))			

			sf = model$coefficient%>%as.numeric
			rf = ((model%>%residuals)^2)%>%mean
			c(sf,rf)			
		}
	
		result = mlapply((1:length(observed)),function(i){
		
		alpha_fit(i)})%>%do.call(rbind,.)
		sfit = result[,1]


		rfit=  result[,2]
		value = mean(rfit^2)
		print(value)
		if(RESBREAK &abs((value_before - value)/value_before) < threshold/10) {print("ResidualBreak");break}

				
		yem = do.call(c,observed)
		
	
		
		xalpha = mlapply(1:length(observed),function(subj_index){
			
			timei1 =timepoints1[[subj_index]]
			timei2 =timepoints2[[subj_index]]
			timei0 =timepoints0[[subj_index]]
			
			xmati  <- lapply(1:length(timei1),function(t){	
			evalb1 = eval.basis(timei1[t], basis1)
			evalb2 = eval.basis(timei2[t], basis2)
			evalb0 = eval.basis(timei0[t], basis0)
			
			
			v = c(evalb0 %x% evalb1 %x% evalb2 )
			return(v)
			})%>%do.call(rbind,.)			
			
			
			xmati*sfit[subj_index]

		})%>%do.call(rbind,.)
		
		
yem = as.matrix(yem)
	library(bigmemory)
	library(bigalgebra)
	A = xalpha
	tA =t(A)
	tyem = t(yem)
	tyem = as.big.matrix(tyem,type="double")
	yem = as.big.matrix(yem,type="double")
	A <- as.big.matrix(A,type="double")
	tA <- as.big.matrix(tA,type="double")	
		
		M1 = as.matrix(tA %*% A )
		M2 = as.matrix(tA %*% yem)		
		
		beta1vec = ginv(M1) %*% M2		
		
		nbasis0=basis0$nbasis
		nbasis1=basis1$nbasis
		nbasis2=basis2$nbasis
		
		
		ind = which.max(abs(beta1vec))		
				

		maxv = 0 
		for(ii in 1:(nbasis1*nbasis2*nbasis0)){
			if (abs((beta1vec_before[ii]-beta1vec[ii])) > maxv){
				maxv = abs((beta1vec_before[ii]-beta1vec[ii]))
			}
		}
		
		
		thresh  = maxv
		beta1vec_before = beta1vec

		for(i in 1:nbasis1){
		for(j in 1:nbasis2){
		for (k in 1:nbasis0){
			beta1[i,j,k] = vectobeta0(beta1vec,i,j,k)
		}
		}
		}
				
		pc_fit = fd2d(beta1, basis1,basis2,basis0)
		pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
		beta1 = pc_fit$betas
		
		thresh =  max(abs(c(beta1_before)-c(beta1)))
		it = it+1
		 {print("Iteration");print(it);print(as.numeric(thresh));print(as.numeric(value))}


	}##while end
	
	
	return(list(beta=beta1,pc_fit = pc_fit,sfit=sfit,thresh = thresh,it=it,value=value))
	
}


if(FALSE){

beta3 = beta1
pc_index=2
observed = observedcenter
result_list[[1]] = res_first
previous_beta[[1]] = res_first$beta
pc_list[[1]] = res_first$pc_fit
betalist =previous_beta

}		

third_FPC_conditional_2d_image  = function(beta3, pc_index, observed,timepoints1,timepoints2,timepoints0,basis1,basis2,basis0, betalist  , threshold=1e-5,minit=1){

library(parallel)
	mlapply<-function(a,f){
	x =mclapply(a,f,mc.cores=60)
	return(x)
	}

start_time <- Sys.time()


thresh = 1
it = 1
value = -1e14

pc_fit = fd2d(beta3, basis1,basis2,basis0)
pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
beta3 = pc_fit$betas
	

pc_fits_previous = lapply(1:length(betalist), function(x){
pc_fit1 = fd2d(betalist[[x]], basis1,basis2,basis0)
pc_fit1$betas = 1/sqrt(norm.fd2d(pc_fit1))*pc_fit1$betas
pc_fit1
})


observed2 =observed
observed2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints1.2 = timepoints1
timepoints1.2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints2.2 = timepoints2
timepoints2.2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints0.2 = timepoints0
timepoints0.2[which(sapply(observed,length)<=pc_index)]=NULL



timei1 =timepoints1.2[[1]]
timei2 =timepoints2.2[[1]]
timei0 =timepoints0.2[[1]]


		
iterationxmatlist_previous=mlapply(1:length(observed),function(subj_index){
	
		iterationxmatcbind = lapply(pc_fits_previous, function(xfit){
		timei1 =timepoints1[[subj_index]]
		timei2 =timepoints2[[subj_index]]
		timei0 =timepoints0[[subj_index]]
			
		iterationxmat = eval.fd2d(timei1,timei2,timei0,xfit)
		return(iterationxmat)
		})%>%do.call(cbind,.)
		
		return(iterationxmatcbind)		
		
})
		




cmat = lapply(pc_fits_previous,constraint.fd2d)%>%do.call(rbind,.)



convcount = 0		
while (thresh >  threshold|it<minit){
	
	
start_time <- Sys.time()


	beta3_before = betatovec(beta3)
	value_before = value
	
	

	iterationxmatlist=mlapply(1:length(observed),function(subj_index){

	timei1 =timepoints1[[subj_index]]
	timei2 =timepoints2[[subj_index]]
	timei0 =timepoints0[[subj_index]]
		
	iterationxmat = eval.fd2d(timei1,timei2,timei0,pc_fit)
	return(iterationxmat)
	})

	
	alpha_fit = function(subj_index){
		
		xmati = iterationxmatlist[[subj_index]]
		xmat_previous = iterationxmatlist_previous[[subj_index]]
								
		model = lm.fit(y=as.matrix(observed2[[subj_index]]), x=cbind(xmat_previous,xmati))
				
		sf = model$coefficient%>%as.numeric
		rf = ((model%>%residuals)^2)%>%mean
		c(sf,rf)
	}

	result= lapply(1:length(observed2),alpha_fit)%>%do.call(rbind,.)
	sfit = result[,1:pc_index]
	rfit=  result[,ncol(result)]
	value = as.numeric(mean(rfit^2))

	
	
	if(RESBREAK & abs((value_before - value)/value_before) < threshold/10) {print('ResidualBreak');break;}
	
	
	N = sapply(observed2,length)%>%sum

	#memory issue
	
	yem = mlapply(1:length(observed2), function(subj_index){
		
	yfits_previous = lapply(1:length(pc_fits_previous), function(pcind){
		sfit[subj_index,pcind]*iterationxmatlist_previous[[subj_index]][,pcind]%>%as.numeric
		})%>%do.call(rbind,.)%>%colSums

	
		
	(observed2[[subj_index]] - yfits_previous)/sqrt(N)
	})%>%do.call(c,.)

	
		
		
		
	xalpha = mlapply(1:length(observed2),function(subj_index){
	
			timei1 =timepoints1[[subj_index]]
			timei2 =timepoints2[[subj_index]]
			timei0 =timepoints0[[subj_index]]
			
			
			xmati  <- lapply(1:length(timei1),function(t){	
			evalb1 = eval.basis(timei1[t], basis1)
			evalb2 = eval.basis(timei2[t], basis2)
			evalb0 = eval.basis(timei0[t], basis0)
			
			
			v = c(evalb0 %x% evalb1 %x% evalb2 )
			return(v)
			})%>%do.call(rbind,.)			
			

			(xmati*sfit[subj_index,ncol(sfit)])/sqrt(N)

	})%>%do.call(rbind,.)
	
	
	A = xalpha
	
	

	tA =t(A)
	tyem = t(yem)
	if(FALSE){
	library(bigmemory)
	library(bigalgebra)
	tyem = as.big.matrix(tyem,type="double")
	A <- as.big.matrix(A,type="double")
	tA <- as.big.matrix(tA,type="double")
	qmat <- 2*as.matrix(tA %*% A)
	pmat = (-2)*as.numeric(as.matrix(tyem%*%A))
	}
	
	qmat <- 2*as.matrix(tA %*% A)
	pmat = (-2)*as.numeric(as.matrix(tyem%*%A))
		
	
	beta3  = lsei::qp(qmat, pmat, cmat, d=rep(0,length(betalist)))	
	max(abs(beta3))
	ind = which.max(abs(beta3))
	
	if (abs(beta3_before[ind]) -abs(beta3[ind]) < 100*threshold){
		if(beta3_before[ind]*beta3[ind] < 0 ){
			beta3= (-1)* beta3
		}
	}
	
	
	threshbefore = thresh
	
	
	nbasis0=basis0$nbasis
	nbasis1=basis1$nbasis
	nbasis2=basis2$nbasis

	maxv = 0 
	for(ii in 1:(nbasis1*nbasis2*nbasis0)){
		if (abs((beta3_before[ii]-beta3[ii])) > maxv){
			maxv = abs((beta3_before[ii]-beta3[ii]))
		}
	}
	thresh  = maxv
	
	if(FALSE){
	if (thresh > threshbefore){convcount = convcount + 1}
	
	if(convcount > 5){
		break;
	}
	}
	
	beta3_before = beta3
	
	beta3 = vectobeta(beta3)
	
	pc_fit = fd2d(beta3, basis1,basis2,basis0)
	pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
	beta3 = pc_fit$betas
    
	it = it+1
	
	{print("Iteration");print(it);print(thresh);print(value)}
	
	

}#while end

return(list(beta=beta3,pc_fit = pc_fit, sfit = sfit, previous_beta =betalist, thresh = thresh,it=it,value=value,gamma=gamma))

}




fpcscore_conditional_2d_image  = function( observed,timepoints1,timepoints2,timepoints0,basis1,basis2,basis0, betalist){

library(parallel)
	mlapply<-function(a,f){
	x =mclapply(a,f,mc.cores=60)
	return(x)
	}

thresh = 1
it = 1
value = -1e14

pc_fits_previous = lapply(1:length(betalist), function(x){
pc_fit1 = fd2d(betalist[[x]], basis1,basis2,basis0)
pc_fit1$betas = 1/sqrt(norm.fd2d(pc_fit1))*pc_fit1$betas
pc_fit1
})


observed2 =observed

		
iterationxmatlist_previous=mlapply(1:length(observed),function(subj_index){
	#print(subj_index)
		iterationxmatcbind = lapply(pc_fits_previous, function(xfit){
		timei1 =timepoints1[[subj_index]]
		timei2 =timepoints2[[subj_index]]
		timei0 =timepoints0[[subj_index]]
			
		iterationxmat = eval.fd2d(timei1,timei2,timei0,xfit)
		return(iterationxmat)
		})%>%do.call(cbind,.)
		
		return(iterationxmatcbind)		
		
})
		


	alpha_fit = function(subj_index){
				
		xmat_previous = iterationxmatlist_previous[[subj_index]]
								
		model = lm.fit(y=as.matrix(observed2[[subj_index]]), x=(xmat_previous))
				
		sf = model$coefficient%>%as.numeric
		rf = ((model%>%residuals)^2)%>%mean
		c(sf,rf)
	}

	result= lapply(1:length(observed2),alpha_fit)%>%do.call(rbind,.)
	sfit = result[,1:length((betalist))]
	rfit=  result[,ncol(result)]
	
	return(sfit)
		
}


pred_soap_2d_fast= function(ylist,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){


timei1 =timepoints1[[1]]
timei2 =timepoints2[[1]]


xmat = lapply(1:length(previous_beta), function(i){
	pc_fit  = pc_list[[i]]
	
	eval.fd2d.image(1:length(timei1),pc_fit)
	
	})%>%do.call(cbind,.)
	
	
cfits = lapply(1:length(ylist), function(x){
	
	index_pc = length(previous_beta)
	
	model = lm.fit(x=as.matrix(xmat[,1:index_pc]),y=as.matrix(ylist[[x]]))
	
	cfit= model%>%coef%>%as.numeric
	cfit
})%>%do.call(rbind,.)

list(sfit = cfits)

}





pred_soap_2d = function(ylist,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){


timei1 =timepoints1[[1]]
timei2 =timepoints2[[1]]

cfits = lapply(1:length(ylist), function(x){
	
	xmat = lapply(1:length(previous_beta), function(i){
		pc_fit  = pc_list[[i]]
		
		eval.fd2d.image(1:length(timei1),pc_fit)
		
	})%>%do.call(cbind,.)
	
	index_pc = length(previous_beta)
	
	cfit = lm(ylist[[x]]~0+xmat[,1:index_pc])%>%coef%>%as.numeric
		
	cfit
})%>%do.call(rbind,.)

yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	
	yfit = lapply(1:length(pc_list),function(index){
		eval.fd2d.image(1:length(timei1),pc_list[[index]])*cfit[index]
	})%>%Reduce('+',.)
	yfit
}
)%>%do.call(cbind,.)


residuals = sapply(1:length(ylist), function(x){
	
	resid = ylist[[x]] - yfits[,x]
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
	
})
sigmahat = mean(residuals[!is.na(residuals)])

list(sigmahat = sigmahat, yfit=  yfits,sfit = cfits)

}

pred_2d<-function(sfits,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){

	cfits = sfits 

	timei1 = timepoints1[[1]]

	cfit = cfits[1,]
	
	yfit = lapply(1:length(pc_list),function(index){
		eval.fd2d.image(1:length(timei1),pc_list[[index]])*cfit[index]
	})%>%Reduce('+',.)		
	
	return (yfit)

}

vectomatrix<-function(vec,picturesize=28){

	firstmat = matrix(0,ncol=picturesize,nrow=picturesize)	
	for (x in 1:picturesize){
	for (y in 1:picturesize){
		
		index= intersect(which(timepoints1[[1]]==x),which(timepoints2[[1]]==y))
		firstmat[x,y] = as.integer(vec[index])
		if (firstmat[x,y]<0){
			firstmat[x,y]=0
		}
		if (firstmat[x,y]>255){
			firstmat[x,y]=255
		}
	}
	}
return (firstmat)
}


pred_soap_2d.nonstandard = function(ylist,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){

timei1 =timepoints1
timei2 =timepoints2

xmatlist  <- lapply(1:length(timei1),function(t){
	
	evalb1 = eval.basis(timei1[t], basis1)
	evalb2 = eval.basis(timei2[t], basis2)
	(t(evalb1) %*% evalb2)				
	}
)



cfits = lapply(1:length(ylist), function(x){
	
	xmat = lapply(1:length(previous_beta), function(i){
		pc_fit  = pc_list[[i]]
		
		eval.fd2d.image.nonstandard2(1:length(timei1),pc_fit,xmatlist)
		
	})%>%do.call(cbind,.)
	
	index_pc = length(previous_beta)
	
	cfit = lm(ylist[[x]]~0+xmat[,1:index_pc])%>%coef%>%as.numeric
		
	cfit
})%>%do.call(rbind,.)

yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	
	yfit = lapply(1:length(pc_list),function(index){
		eval.fd2d.image.nonstandard2(1:length(timei1),pc_list[[index]],xmatlist)*cfit[index]
	})%>%Reduce('+',.)		
	yfit
}
)%>%do.call(cbind,.)


residuals = sapply(1:length(ylist), function(x){
	
	resid = ylist[[x]] - yfits[,x]
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
	
})
sigmahat = mean(residuals[!is.na(residuals)])

list(sigmahat = sigmahat, yfit=  yfits,sfit = cfits)

}


