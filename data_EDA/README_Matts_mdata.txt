
M- matrix of age by year where all values are 0.2
weight-taken from stock_wt.dat, weights for ages 2-14, data for 1959-1982 are based on the median weight for each age from the available data
mat- taken from mat.txt, subset for ages 2-14
midy_weight- taken from midy_wt.dat, weights for ages 2-14, data for 1959-1982 are based on the median weight for each age from the available data
comm_wt- weights for age 3-14. Comes from a csv file that I created (comm_wt_mod) that was based on comm_wt.txt where 2015 and 2016 weights are calculated based on the mean of the past 5 years
C- catch for ages 2-14+, 14+ is just a sum of ages 14, 15, 16 (NAs treated as 0's for the sum), all zeros treated instead as low values (0.0005). This was in Noel's code
log_C- just the log of C
landings- Noel's code had landings based on catch*midywt. I just changed this so that landings come from landings.dat instead. Additionally, landings.dat had two values for 2016, I summed these values. Will need to ask Noel why there are 2 values 
log_landings- just the log of landings
index- 1983-2016 for ages 2-14. Values based on a new index calculated by multiplying OFF values by the mean ratio between IO and OFF surveys. I believe this is what was used in the last stock assessment. Can easily be switched to use OFF or IO alone by changing indices in mdata to the corresponding survey vector
