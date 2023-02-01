

#for (i in 1 : n){ 
#    
#    for(j in 1 : m){ # this for-loop should be paralellized
#        
#        # some code that can be done in parallel
#        
#    }
#    
#    # some more code
#    
#}


#
# Reproducible example

n <- 100
m <- 1000
# in my code, n=100,000 and m=40,000
output <- numeric(n) # vector of output
tmp    <- numeric(m) # temporary object

for (i in 1 : (n-1)){ 
    
    for(j in 1 : m){ # this for-loop should be paralellized
        
        tmp[j] <- rnorm(1, mean = output[i], sd = 1) # can be done in parallel
        
    }
    
    output[i+1] <- mean(tmp)
    
}
