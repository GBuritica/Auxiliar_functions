# Auxiliar_functions
- Read 1st
 ## Marginal transformation to Fréchet domain of attraction. 
 ## Ranks transform  : Given a trajectory of observations:
 ##                      ranks              -> returns the vector of ranks
 ##                      ranktranform  -> returns the vector with rank transform applied

Examples of implementation: 

 sample  <- arima.sim(n = 5000, list(ar=0.5, ma=0.3), rand.gen=function(n) rt(n,df=1))
 plot.ts(abs(sample))
 plot.ts(ranktransform(sample2))
 ranktransform(c(NA,NA,3,5,5))


General comments: 
You might want to handle repeated values and NA's before applying the functions.

Comments to each function: 

ranks: NA values are transformed to 0 and then sorted in order of appearence.
ranktranform: NA and 0's set to 0. 
