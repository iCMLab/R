
###########################################################################################
#               Msc. PELAYO G. DE LENA RODR√çGUEZ                                          #
#               R version 3.5.3 (2019-03-11) 
#                 RStudio Version 1.2.1335
#         custom-script; cp all plots in current session R
#     Instituto de Estudios Celulares y Moleculares - ICM                                #
#                     MADRID, 16/05/2019                                                 #
###########################################################################################

# # COPIAR PLOTS EN ACTUAL SESIÛN



plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="./")


