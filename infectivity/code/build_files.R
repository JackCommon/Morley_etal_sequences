phage <- seq(1,12) %>% 
  rep(3)
host <- phage

replicate <- rep("2.7", 36)
timepoint <- c( rep("t1", 12),
                rep("t4", 12),
                rep("t9", 12))

infectivity <- rep(0,36)
resistance  <- rep(0,36)

seven <- data_frame(replicate, timepoint, phage, host, infectivity, resistance)

write.csv(file="./time_shift/summary_data/infectivity_data_parts/2.7.csv", seven,
          row.names = F)
