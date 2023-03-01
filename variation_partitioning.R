library(vegan)
library(adespatial)


space_fw <- forward.sel(y, space,nperm=999, alpha = 1,R2more= 0.01, verbose=F)   
selected.space <- space[, space_fw$order]

soil_fw <-forward.sel(y, soil,nperm=999, alpha = 1,R2more= 0.01, verbose=F) 
selected.soil <- soil[, soil_fw$order]

vegetation_fw <-forward.sel(y, vegetation,nperm=999, alpha = 1,R2more= 0.01,verbose=F)    
selected.veg <- vegetation[, vegetation_fw$order]

temporal_fw <-forward.sel(y, temporal,nperm=999, alpha = 1,R2more=0.01, verbose=F)   
selected.temp <- temporal[, temporal_fw$order]

vartpart_a.out <- varpart( y  ,selected.veg,selected.soil,selected.space,selected.temp)

#########

floristic_fw=forward.sel(y, floristic,nperm=999, alpha = 1,R2more= 0.01, verbose=F) #  
selected.floristic = floristic[, floristic_fw$order]

host_phylogeny_fw=forward.sel(y, host_phylogeny,nperm=999, alpha = 1,R2more= 0.01,verbose=F) #  
selected.host_phylogeny= host_phylogeny[, host_phylogeny_fw$order]

host_sp_fw=forward.sel(y, host_sp,nperm=999, alpha = 1,R2more=0.01, verbose=F) 
selected.host_sp= host_sp[, host_sp_fw$order]

others_fw=forward.sel(y, others,nperm=999, alpha = 1,R2more= 0.01, verbose=F)  
selected.others= others[, others_fw$order]
 
vartpart_b.out =varpart(y, selected.host_sp, selected.floristic, selected.host_phylogeny,selected.others)
 
