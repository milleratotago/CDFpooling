# Make an example file for testing
# This example has 2 conditions with RTs in cond2 a bit slower than those in cond1.
# It is not intended to generate realistic RT data but only to generate
# two-condition example data in the correct format for the CDFpooling routines.

library(gamlss)

file_name <- "CDFpooling_example.csv"
n_conditions <- 2
n_trials_per_sub_cond <- 200
n_subs <- 10

n_trials_ttl <- n_subs * n_conditions * n_trials_per_sub_cond

mu <- 600
sigma <- 60
tau <- 100

rts <- round(rexGAUS(n_trials_ttl,mu,sigma,tau))

sub <- sample(1:n_subs,n_trials_ttl,replace=T)
cond <- sample(1:n_conditions,n_trials_ttl,replace=T)
rt_df <- data.frame(sub=sub,cond=cond,rt=rts,row.names=NULL)
rt_df <- rt_df[order(sub),]
cond2 <- cond == 2
rt_df$rt[cond2] <- round( 1.1*rt_df$rt[cond2] + 20 )
write.csv(rt_df, file = file_name, row.names = FALSE)
