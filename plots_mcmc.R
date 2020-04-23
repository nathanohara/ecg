library(ggplot2)
library(dplyr)
library(tidyverse)

pred_mean <- unlist(map(MEAN, mean))

MEAN_burn <- map(seq(MEAN), function(x) MEAN[[x]][1000:1707])

pred_mean <- unlist(map(MEAN_burn, mean))


pred_df = tibble(
  time = t,
  pred_rate = pred_mean,
  emp_rate = y
)

pred_df %>% 
  ggplot() +
  geom_line(aes(x = time, y = emp_rate), col = "blue") +
  geom_line(aes(x = time, y = pred_rate), col = "orange")

samples %>%
  ggplot() +
  geom_line(aes(x = time, y = signal), col = "blue") 

par(mfrow = c(2, 1))
for(j in sample(12, 4)){
  plot(ALPHA[[j]], type = 'l', ylab = paste("Alpha", j), xlab = "Iter", col = "purple",
       main = paste("Acceptance ratio =", round(length(unique(ALPHA[[j]])) / S, 3)))
  acf(ALPHA[[j]], ylab = paste("Alpha", j), main = paste("Series Alpha", j))
}


plot(DELTA[[3]], type = 'l', ylab = paste("Delta", 3), xlab = "Iter", col = "purple",
     main = paste("Acceptance ratio =", round(length(unique(DELTA[[3]])) / S, 3)))
acf(DELTA[[3]], ylab = paste("Delta", 3), main = paste("Series Delta", 3))

for(j in sample(17, 2)){
  plot(DELTA[[j]], type = 'l', ylab = paste("Delta", j), xlab = "Iter", col = "purple",
       main = paste("Acceptance ratio =", round(length(unique(DELTA[[j]])) / S, 3)))
  acf(DELTA[[j]], ylab = paste("Delta", j), main = paste("Series Delta", j))
}

par(mfrow = c(1,2))
plot(DC, type = 'l', ylab = "DC", xlab = "Iter", col = "purple",
       main = paste("Acceptance ratio =", round(length(unique(DC)) / S, 3)))
acf(DC, ylab = "DC", main = "Series DC")

par(mfrow = c(1,2))
plot(TAO, type = 'l', ylab = "Tao", xlab = "Iter", col = "purple",
    TAOin = paste("Acceptance ratio =", round(length(unique(TAO)) / S, 3)))
acf(TAO, ylab = "Tao", main = "Series Tao")



save.image(file = "Trial_2.RData")
dir()
load("Trial_2.RData")
