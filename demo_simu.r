library(dplyr)
library(fda)
source('funs_mdfpca.r')
load('brain_fpc.RData')

load('modelfp')

basetimepoints1 = c()
basetimepoints2 = c()
for (i in 1:30){
for (j in 1:30){
	basetimepoints1 = c(basetimepoints1,i)
	basetimepoints2 = c(basetimepoints2,j)	
}
}

basetimepoints0 = c(0,5,10,15,20,25,30)


gentimepoints<-function(timepoints0unique){

timepoints0=c()
timepoints1=c()
timepoints2=c()
for(i in 1:length(timepoints0unique)){
timepoints0 = c(timepoints0,rep(timepoints0unique[i],30*30))
timepoints1 = c(timepoints1,basetimepoints1)
timepoints2  = c(timepoints2,basetimepoints2)
}

return(list(timepoints0=timepoints0,timepoints1=timepoints1,timepoints2=timepoints2))
}



result_list_global<<-result_list


allresultlist=list()

for(sim in 1:50){

varlistmat =lapply(1:100,function(i){
var1= rnorm(1,0,(4000))
var2= rnorm(1,0,(1000))
var3= rnorm(1,0,(500))
return(c(var1,var2,var3))
})%>%do.call(rbind,.)



olist = lapply(1:100,function(ss){

timepoints0unique = basetimepoints0
temp = gentimepoints(timepoints0unique)
timepoints0 = temp$timepoints0
timepoints1 = temp$timepoints1
timepoints2 = temp$timepoints2

numofsimupc=2

values_global=list()
for(pcindex in 1:numofsimupc){
	valuev = eval.image.process(pc_list[[pcindex]],timepoints0,timepoints1,timepoints2)		
	values_global[[pcindex]]=valuev
}

outputv=c()
varlist = varlistmat[ss,]


error = rnorm(length(timepoints1),0,0.01)
for (index in 1:length(timepoints1)){
	
	e = error[index]
	tempvalue = 0
	
	for (pcindex in 1:numofsimupc){
		tempvalue = tempvalue + values_global[[pcindex]][index]*varlist[pcindex]+e
	}
	
	outputv = c(outputv,tempvalue)
}

return(c(outputv,timepoints0,timepoints1,timepoints2))

})


observedlist = lapply(olist,function(item){
k = length(item)/4
item[(1):(k)]
})

timepoints0list = lapply(olist,function(item){
k = length(item)/4
item[(k+1):(2*k)]
})

timepoints1list = lapply(olist,function(item){
k = length(item)/4
item[(2*k+1):(3*k)]
})

timepoints2list = lapply(olist,function(item){
k = length(item)/4
item[(3*k+1):(4*k)]
})









maxpixelsize = 30
nbasis=12
library(fda)
nbasis0 = 4

basis1=create.bspline.basis(rangeval=c(0,maxpixelsize+1),nbasis=nbasis,norder=4)
basis2=create.bspline.basis(rangeval=c(0,maxpixelsize+1),nbasis=nbasis,norder=4)
basis0=create.bspline.basis(rangeval=result_list_global[[1]]$pc_fit$basis0$rangeval,nbasis=nbasis0,norder=4)
nbasis1=basis1$nbasis
nbasis2=basis2$nbasis
nbasis0=basis0$nbasis




beta1 = array(1,c(nbasis,nbasis,nbasis0))
for(k in 1:nbasis0){
for(i in 1:nbasis){
	beta1[i,,k]=seq(0.001,0.012,by=0.001)
}
}


start.time = proc.time()

previous_beta = list()
pc_list = list()
result_list=list()

numIter = 2
for(i in 1:numIter){

	if(i==1){
		res_first = first_FPC_2d_image(beta1,observedlist,timepoints1list,timepoints2list,timepoints0list,basis1,basis2,basis0,threshold=1e-4,minit=3)
		result_list[[1]] = res_first
		previous_beta[[1]] = res_first$beta
		pc_list[[1]] = res_first$pc_fit
		
	}else{
		
		res_second = third_FPC_conditional_2d_image(beta1, pc_index=i, observedlist,timepoints1list,timepoints2list,timepoints0list,basis1,basis2,basis0, betalist =previous_beta , threshold=1e-4,minit=1)
		result_list[[i]] = res_second
		previous_beta[[i]] = res_second$beta
		pc_list[[i]] = res_second$pc_fit
		

	}

}
proc.time() - start.time







###IMSE
outputv2=c()
for (pcindex in 1:numofsimupc){

	xy = 2*inprod.fd2dx(result_list[[pcindex]]$pc_fit,result_list_global[[pcindex]]$pc_fit)
	x2=inprod.fd2d(result_list[[pcindex]]$pc_fit,result_list[[pcindex]]$pc_fit)
	y2=inprod.fd2d(result_list_global[[pcindex]]$pc_fit,result_list_global[[pcindex]]$pc_fit)
		
	outputv2 =c(outputv2,min(x2+y2-xy,x2+y2+xy))
}

imse = outputv2
fpcscore = result_list[[numofsimupc]]$sfit

for(i in 1:ncol(fpcscore)){
if (mean(abs(fpcscore[,i] - varlistmat[,i])) > mean(abs((-1)*fpcscore[,i] - varlistmat[,i]))){
fpcscore[,i] = fpcscore[,i] * (-1)
}
}

reladiff_fpcscore = sapply(1:ncol(fpcscore),function(i){
mean(abs(fpcscore[,i] - varlistmat[,i]))/sd(varlistmat[,i])
})

o = list(imse = imse, fpcscore=fpcscore, truefpcscore = varlistmat,result_list = result_list)

allresultlist[[sim]] = o 
}

