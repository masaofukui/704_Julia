packages_vector <-c("ggplot2","OECD","data.table","dplyr","tidyr",
                    "countrycode","zoo","haven","readstata13","seasonal",
                    "stringr","nleqslv")
#install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE) # the "lapply" function means "apply this function to the elements of this list or more restricted data type"
rm(list = ls())
setwd("/Users/fukui/Dropbox (Personal)/teaching/704B/stata")
df <- fread("./data/fred_unrate.csv")
fs_solve <- function(x,f_bias,s_bias) {
  y <- numeric(2)
  y[1] <- x[1]*(1- exp(-x[1]-x[2]))/(x[1] + x[2]) - s_bias
  y[2] <- x[2]*(1- exp(-x[1]-x[2]))/(x[1] + x[2]) - f_bias
  y
}
fs <- c(0.2,0.01)
fs_data_in <- function(fs){
  if (!is.na(fs[1]) ){
    xstart = c(0.02,0.1)
    result <- nleqslv(xstart, fs_solve,jac = NULL,fs[1],fs[2])
    fs_sol <- result$x
  }else{
    fs_sol <- c(NA,NA)
  }
  fs_sol
}

a <- df[,c("job_finding_bias","separation_bias")]
a <- as.list(data.frame(t(a)))


result <- lapply(a,fs_data_in)
b <- as.data.frame(do.call(rbind,result))

df_new <- cbind(df,b)
df_new <- df_new %>% dplyr:: rename(separation = V1, job_finding = V2)

fwrite(df_new,"./data/fred_unrate_bias_corrected.csv")
