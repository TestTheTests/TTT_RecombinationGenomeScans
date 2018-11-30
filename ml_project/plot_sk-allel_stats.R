####################################################################################################################
#
# File    : plot_sk-allel_stats.R
# History : 11/20/2018 - Created by Kevin Freeman (KF)
#
#################################################################################################################
#
# Plots of all the stats calculated in sk-allel, used as a sanity-check to make sure they are being calculated
# correctly. The plots include red lines indicating where the inversion is and some of them have blue lines
# for the variable recombination region and purple lines for the sweep
#
#################################################################################################################



resultsDir <- "/media/kevin/TOSHIBA_EXT/TTT_RecombinationGenomeScans/results_final"
sim   <- read.csv(paste0(resultsDir, "/10900_Invers_ScanResults_with_sk-allel.txt"), header = TRUE, sep = " ")

inversion <-c(320000, 330000)
sweep     <- 175000
varRC     <- c(400000, 450000)

plot(sim$pos, sim$tajD_sk.allel_v1.1.10, main = "Tajima's D")
abline(v = inversion, col = "red")
plot(sim$pos, sim$pi_sk.allel_v1.1.10, main = "Pi (sequence diversity)")
abline(v = inversion, col = "red")
plot(sim$pos, sim$thetaW_sk.allel_v1.1.10, main = "Watterson's Theta")
abline(v = inversion, col = "red")
plot(sim$pos, sim$H12_sk.allel_v1.1.10, main = "H12")
abline(v = sweep, col = "purple")
abline(v = varRC, col = "blue")
abline(v = inversion, col = "red")
plot(sim$pos, sim$ihs_sk.allel_v1.1.10, main = "iHS (sk-allel)")
abline(v = inversion, col = "red")
abline(v = varRC, col = "blue")
plot(sim$pos, sim$nsl_sk.allel_v1.1.10, main = "NSL")
abline(v = inversion, col = "red")
abline(v = varRC, col = "blue")
plot(sim$rehh_2.0.2_ALL_iHS, sim$ihs_sk.allel_v1.1.10, main = "iHS (rehh) vs iHS (sk-allel)",
     xlim = c(-4,4), ylim = c(-4,4))
plot(sim$pos, sim$rehh_2.0.2_ALL_iHS, main = "iHS (rehh) vs iHS(sk-allel)")
points(sim$pos, sim$ihs_sk.allel_v1.1.10, col = "red")
