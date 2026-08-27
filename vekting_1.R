#renv::init()
#renv::install("sampling")
#renv::install("statisticsnorway/ReGenesees")
library(sampling)
library(ReGenesees)
load("populasjon.RData")
head(d)

#------------------------------------------------------
# Stratifisert tilfeldig trekking etter kjonn         
# krysset med aldersgruppene 16-24, 25-44, 45-66, 67+ 
# Total utvalgsstørrelse på 1000.                     
# Samme trekkeandel i alle strata.                    
#------------------------------------------------------
# Trekkerammen
U <- d[d$alder >= 16, ]

# Utvalgsstørrelsen i hvert stratum
n <- round(table(U$kjonn, U$ald_grp) * 1000 / nrow(U))
n
U <- U[order(U$kjonn, U$ald_grp), ]

# Trekker utvalg
s <- strata(data = U,
            stratanames = c("kjonn", "ald_grp"),
            size = c(67, 169, 173, 102,
                     64, 162, 165, 98),
            method = "srswor")

# Legger på alle variabler fra trekkerammen
utvalg <- getdata(U, s)

# Designvekten
utvalg$d_vekt <- 1 / utvalg$Prob
head(utvalg)

#----------------------------------------------
# Vi leser inn en tidligere trukket utvalgsfil 
# slik at alle jobber med samme data.          
#----------------------------------------------
load("utvalg.RData")
head(utvalg)

# Estimerer vektet gjennomsnitt og sammenligner med snitt i populasjon
weighted.mean(utvalg$y1, w = utvalg$d_vekt)
mean(U$y1)

weighted.mean(utvalg$y2, w = utvalg$d_vekt)
mean(U$y2)

weighted.mean(utvalg$y3, w = utvalg$d_vekt)
mean(U$y3)

#----------------------------------------------
# Med frafall: Vi leser inn en tidligere laget 
# nettofil, slik at alle jobber med samme data.                   
#----------------------------------------------
load("netto.RData")
head(netto)

#------------------------------------------------------------
# Metode 1: Vi bare justerer designvektene basert på antagelse 
# om helt tilfeldig frafall "Missing Complete at Random" (MCAR)
#------------------------------------------------------------
netto$mcar_vekt <- netto$d_vekt * nrow(utvalg) / nrow(netto)
#
weighted.mean(netto$y1, w = netto$mcar_vekt) 
mean(U$y1)
#
weighted.mean(netto$y2, w = netto$mcar_vekt) 
mean(U$y2)
#
weighted.mean(netto$y3, w = netto$mcar_vekt) 
mean(U$y3)


