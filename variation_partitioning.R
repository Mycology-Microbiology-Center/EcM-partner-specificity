library(vegan)
library(adespatial)

 

space_fw=forward.sel(y, x,nperm=999, alpha = 1,R2more= 0.01, verbose=F)    
selected.space= x[, space_fw$order ]

soil_fw=forward.sel(y, x,nperm=999, alpha = 1,R2more= 0.01, verbose=F)   
selected.soil= x[, soil_fw$order ]

vegetation_fw=forward.sel(y, x,nperm=999, alpha = 1,R2more= 0.01,verbose=F)   
selected.veg= x[, vegetation_fw$order ]


temporal_fw=forward.sel(y, x,nperm=999, alpha = 1,R2more=0.01, verbose=F) 
selected.temp= x[, temporal_fw$order ]

vartpart.out <- varpart( y  ,selected.veg,selected.soil,selected.space,selected.temp)
