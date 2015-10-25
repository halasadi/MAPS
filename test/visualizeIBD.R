four_cm = read.table("data/popres/popres_lowerBnd_4_upperBnd_Inf-EEMS2-test-sim/rdistJtDobsJ.txt")
sam <- sammon(max(four_cm)-as.matrix(four_cm))
plot(sam$points, type = "n")

rownames(four_cm) = c("France", "Ireland", "Italy", "Portugal", "Spain", "United Kingdom")
colnames(four_cm) = row.names

text(sam$points, labels = countries_to_keep)
