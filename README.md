# nonParQuntileCausality

Nonparametric Quantile Causality Test  

lrq.causality.test(x, y, type=c("mean","variance"), q=NULL, hm=NULL)  

x = numeric vector, cause (independent) variable   
y = numeric vector, dependent variable  
q = vector of quatiles, default = 0.01,....,0.99   
type = "mean" for causality in mean test,  
       "variance" for causality in variance test   
hm = numeric scalar, bandwidth,  
     if NULL  optimal bandwith of Yu and Jones (1998) is used  

code only considers first lags of x and y. 

Returns:  
   stat = vector t statistics for causality at each quantile   
   q = vector of quantiles 

Mehmet Balcilar, 2014-7-4  

