library(multimode)
library(mixtools)

args <- commandArgs(trailingOnly = TRUE)
export_path <- args[1]
file_path <- args[2]
row_index <- as.numeric(args[3])

df = read.csv("../../pipeline_testing/distance_and_ordering/SMC01_N/cell_distances.csv", row.names=1)
df = as.numeric(df[row_index,])

# #Testing multimodality (Null hypothesis is - unimodality, alternate hypothesis - multimodality)
#print(modetest(df))

# 1. Calculate the optimal bandwidth for your data
bandwidth <- bw.nrd0(df)

# 2. Calculate the number of modes directly
num_modes <- nmodes(df, bw = bandwidth)

# #finding modes (Maximum frequency) and antimodes (least frequent value between two modes)
result <- locmodes(df,mod0=num_modes,display=TRUE,posLegend="topright")
print(result$locations)

## Fitting the data with more than one normal distributions
fit <- normalmixEM(df,k=num_modes)
mu <- fit$mu
sigma <- fit$sigma
lambda <- fit$lambda
fit_result <- as.data.frame(cbind(mu,sigma,lambda))

output_file_name <- paste(export_path,"norm_dist_param_R.csv",sep="/")

write.csv(fit_result, file = output_file_name, row.names = FALSE)