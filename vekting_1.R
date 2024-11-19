#renv::init()
#renv::install("sampling")
#renv::install("statisticsnorway/ReGenesees")
library(sampling)
library(ReGenesees)
load("trekkeramme.RData")
head(d)

#######################################################
# Stratifisert tilfeldig trekking etter kjonn         #
# krysset med aldersgruppene 16-24, 25-44, 45-66, 67+ #
# Total utvalgsstørrelse på 1000.                     #
# Samme trekkeandel i alle strata.                    #
#######################################################
# Trekkerammen
U <- d[d$alder >= 16, ]
# Aldersinndelingen
U$ald_grp <- cut(U$alder, breaks = c(15, 24, 44, 66, 100))
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

################################################
# Vi leser inn en tidligere trukket utvalgsfil #
# slik at alle jobber med samme data.          #
################################################
load("utvalg.RData")
head(utvalg)
# Estimerer populasjonstotal (i millioner) med Horvitz-Thomsen-estimatoren
HTestimator(utvalg$y1, utvalg$Prob) / 1e+6
# Den faktiske totalen er:
sum(U$y1) / 1e+6
#
HTestimator(utvalg$y2, utvalg$Prob) / 1e+6
sum(U$y2) / 1e+6

#########################################################
# Med frafall: Vi leser inn en tidligere laget nettofil #
# slik at alle jobber med samme data.                   #
#########################################################
load("netto.RData")
head(netto)
# Vi justerer designvektene basert på antagelse om helt tilfeldig frafall
# "Missing Complete at Random" (MCAR)
netto$mcar_vekt <- netto$d_vekt * nrow(utvalg) / nrow(netto)
#
crossprod(netto$mcar_vekt, netto$y1) / 1e+6 
sum(U$y1) / 1e+6
#
crossprod(netto$mcar_vekt, netto$y2) / 1e+6
sum(U$y2) / 1e+6
#
row.names(netto) <- NULL
save(netto, file = "netto.RData" )

##################################
# Etterstratifisering mot region #
##################################
t_U <- table(U$region)
t_n <- table(netto$region)
w <- t_U / t_n
w
netto$es_vekt <- w[netto$region]
#
# Estimering av populasjonstotaler (vektet sum)
crossprod(netto$es_vekt, netto$y1) / 1e+6
sum(U$y1) / 1e+6
#
crossprod(netto$es_vekt, netto$y2) / 1e+6
sum(U$y2) / 1e+6

################################################
# Etterstratifisering mot kjonn x aldersgruppe #
################################################
t_U <- table(U$kjonn, U$ald_grp)
t_n <- table(netto$kjonn, netto$ald_grp)
w <- t_U / t_n
w
# Den etterstratifiserte vekten
netto$es_vekt <- w[cbind(netto$kjonn, netto$ald_grp)]
#
# Estimering av populasjonstotaler (vektet sum)
crossprod(netto$es_vekt, netto$y1) / 1e+6
sum(U$y1) / 1e+6
#
crossprod(netto$es_vekt, netto$y2)
sum(U$y2) / 1e+6

############################################################################
# Kalibrering med ReGenesees. Kalibreringsmodell: region + kjonn * ald_grp #
############################################################################
des <- e.svydesign(data = netto, ids = ~pnr, weight = ~d_vekt)
tot_temp <- pop.template(data = des,
                         calmodel = ~region + kjonn * ald_grp - 1)
tot_temp
tot_pop <- fill.template(universe = U,
                         temp = tot_temp)
tot_pop
descal <- e.calibrate(design = des,
                      df.population = tot_pop,
                      calfun = "linear")
netto$kal_vekt <- weights(descal)
head(netto)
#
summary(netto$d_vekt / mean(netto$d_vekt) * mean(netto$kal_vekt) / netto$kal_vekt)
##################################
svystatTM(descal, 
          y = ~y1, 
          vartype = "se") / 1e+6

sum(U$y1) / 1e+6
######################################
svystatTM(descal, 
          y = ~y2, 
          vartype = "se") / 1e+6

sum(U$y2) / 1e+6
######################################
svystatTM(descal, 
          y = ~y1, 
          by = ~region, 
          vartype = "se")

tapply(U$y1, U$reg, sum)
#######################################
svystatTM(descal, 
          y = ~y1, 
          by = ~kjonn * ald_grp, 
          vartype = "se")

tapply(U$y1, list(U$kjonn, U$ald_grp), sum)
#############################################
svystatTM(descal, 
          y = ~y1, 
          by = ~utd, 
          estimator = "Mean", 
          vartype = "se")

tapply(U$y1, U$utd, mean)
############################################
svystatL(descal, 
         expression(y1 / y2), 
         vartype = "se")

sum(U$y1) / sum(U$y2)
############################################
svystatL(descal, 
         expression(y1 / y2), 
         by = ~region, 
         vartype = "se", 
         conf.int = TRUE)

tapply(U$y1, U$region, sum) / tapply(U$y2, U$region, sum)
##########################################################

  
