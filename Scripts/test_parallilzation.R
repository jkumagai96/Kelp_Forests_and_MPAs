library(doParallel)
library(foreach)

print("packages loaded")

# Start the cluster for parallel processing 
n_cores <- 5
my.cluster <- makePSOCKcluster(n_cores)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParWorkers()

a <- c(1:1000)
b <- c(5)
print('b')

start <- Sys.time()
bootstrap_list <- foreach (j = 1:1000) %dopar% {

values <- a*b
values

}
end <- Sys.time()

start - end

print('end of script')
