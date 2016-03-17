#Intsaller la librairie testthat dans le R
#A lancer depuis le R avec la commande testthat::test_dir("chemin_des_sources_du_package/tests")
library(testthat)
library(VAM)

test_check("VAM")
