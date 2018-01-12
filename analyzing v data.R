# Reading data ------
datasets <- vector(5, mode = "list")
for(i in 1:5) {
  dat <- read.csv(paste("data/data_1_", data_num, ".csv", sep = ""))
  dat$dataset <- i
  datasets[[i]] <- dat
}
