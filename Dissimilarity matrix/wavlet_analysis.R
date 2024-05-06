library(WaveletComp)
# Specify the path to the CSV file

save_path <- "../output/"

# Load the CSV file into a data frame
df<- read.csv("../output/data_ww_geodata_CA_cdph_scan.csv") 

ts_data <- ts(df$SC2_N_gc_g_dry_weight, start = 1, end = length(df$SARS_cov2), frequency = 1)

#x = periodic.series(start.period = 50, length = 1000)
#x = x + 0.2*rnorm(1000) # add some noise

#my.data <- data.frame(x = x)
#my.w <- analyze.wavelet(my.data, "x",
#                        loess.span = 0,
#                        dt = 1, dj = 1/250,
#                        lowerPeriod = 16,
#                        upperPeriod = 128,
#                        make.pval = TRUE, n.sim = 10)
#wt.image(my.w, color.key = "quantile", n.levels = 250,
#         legend.params = list(lab = "wavelet power levels", mar = 4.7))

unique_plants <- unique(df$Plant)

i=3
j=2
plant1 <- unique_plants[i]
plant2 <- unique_plants[j]
  
  # Extract time series data for the two plants
ts1_ <- df[df$Plant == plant1, "SC2_N_gc_g_dry_weight"]

ts1<- ts1_[11:length(ts1_)]
ts2 <- df[df$Plant == plant2,"SC2_N_gc_g_dry_weight"]

my.data <- data.frame(x = ts1, y = ts2) 
my.wc <- analyze.coherency(my.data, my.pair = c("x","y"))
wc.image(my.wc,n.levels=250, legend.params=list(lab="cross-waveletpowerlevels"), timelab="",periodlab="period(days)")
