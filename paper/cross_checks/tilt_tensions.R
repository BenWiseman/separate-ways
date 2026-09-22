# tilt_tensions.R -- the Turok-Boyle tilt n_s = 1 - 7 alpha_3(M_P)/pi against 2026 CMB determinations
# (Paper 2 sec 5.2 table). Base R only. Inputs: the table's quoted (n_s, sigma) per dataset, verbatim
# from SUPPLEMENT_v3.md sec S6 (which quotes P8_REPORT.md); the prediction from alpha_3(M_P) = 0.0189 [41].
source("helpers.R")
alpha3 <- 0.0189; ns_pred <- 1 - 7*alpha3/pi
cat(sprintf("n_s(pred) = 1 - 7*%.4f/pi = %.6f (paper: 0.957888)\n", alpha3, ns_pred))
tab <- data.frame(
  dataset = c("Planck 2018 TT,TE,EE+lowE+lensing","ACT DR6 alone","SPT-3G D1 alone","P-ACT","CMB-SPA","P-ACT-LB","P-ACT-LB2","CMB-SPA + DESI DR2"),
  ns    = c(0.9649,0.9666,0.9510,0.9709,0.9679,0.9743,0.9752,0.9726),
  sigma = c(0.0042,0.0077,0.0110,0.0038,0.0033,0.0034,0.0030,0.0028),
  paper = c(1.67,1.13,-0.63,3.42,3.03,4.83,5.77,5.25))
tab$tension <- (tab$ns - ns_pred)/tab$sigma
print(tab, digits=4)
cat("\nchecks:\n")
report("n_s prediction", expected = 0.957888, reproduced = ns_pred, tol = 1e-5, mode = "rel")
for (i in seq_len(nrow(tab)))
  report(sprintf("tilt tension, %s [sigma]", tab$dataset[i]), expected = tab$paper[i], reproduced = tab$tension[i], tol = 0.02, mode = "rel",
         note = "paper sec 5.2 table; inputs quoted there")
