setwd("/ssb/bruker/mld/kurs_utvalg")
renv::init()
renv::install("sampling")
renv::install("statisticsnorway/ReGenesees")
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
U <- d[d$alder >= 16, ]
U$ald_grp <- cut(U$alder, breaks = c(15, 24, 44, 66, 100))
# Utvalgsstørrelse i hvert stratum
n <- round(table(U$kjonn, U$ald_grp) * 1000 / nrow(U))
n
U <- U[order(U$kjonn, U$ald_grp), ]
s <- strata(data = U,
            stratanames = c("kjonn", "ald_grp"),
            size = c(n[1,], n[2, ]),
            method = "srswor")
utvalg <- getdata(U, s)
utvalg$d_vekt <- 1 / utvalg$Prob
head(utvalg)
#
HTestimator(utvalg$y1, utvalg$Prob)
sum(U$y1)
#
HTestimator(utvalg$y2, utvalg$Prob)
sum(U$y2)

###############
# Med frafall #
###############
# Responssannsynlighet
utvalg$rs <- ifelse(utvalg$utd == "4", 0.60,
                    ifelse(utvalg$utd == "3", 0.50,
                           ifelse(utvalg$utd == "2", 0.40, 0.30)))
# Nettoutvalget
netto <- utvalg[utvalg$rs > runif(nrow(utvalg)), ]

table(netto$kjonn, netto$ald_grp)

# Ny vekt basert på antagelse om helt tilfeldig frafall (MCAR)
netto$mcar_vekt <- netto$d_vekt * nrow(utvalg) / nrow(netto)
#
crossprod(netto$mcar_vekt, netto$y1)
sum(U$y1)
#
crossprod(netto$mcar_vekt, netto$y2)
sum(U$y2)

#############################################
# Etterstratifisering: kjonn x aldersgruppe #
#############################################
t_U <- table(U$kjonn, U$ald_grp)
t_n <- table(netto$kjonn, netto$ald_grp)
w <- t_U / t_n
w
# Etterstratifisert vekt
netto$ps_vekt <- w[cbind(netto$kjonn, netto$ald_grp)]
#
crossprod(netto$ps_vekt, netto$y1)
sum(U$y1)
#
crossprod(netto$ps_vekt, netto$y2)
sum(U$y2)

##########################################
# Etterstratifisering: kjonn x utdanning #
##########################################
t_U <- table(U$kjonn, U$utd)
t_n <- table(netto$kjonn, netto$utd)
w <- t_U / t_n
w
netto$ps_vekt <- w[cbind(netto$kjonn, netto$utd)]
#
crossprod(netto$ps_vekt, netto$y1)
sum(U$y1)
#
crossprod(netto$ps_vekt, netto$y2)
sum(U$y2)

###############################
# Etterstratifisering: region #
###############################
t_U <- table(U$region)
t_n <- table(netto$region)
w <- t_U / t_n
w
netto$ps_vekt <- w[netto$region]
#
crossprod(netto$ps_vekt, netto$y1)
sum(U$y1)
#
crossprod(netto$ps_vekt, netto$y2)
sum(U$y2)

##############################################################
# Kalibrering med ReGenesees: region + utd + kjonn * ald_grp #
##############################################################
des <- e.svydesign(data = netto, ids = ~pnr, weight = ~Prob)
tot_temp <- pop.template(data = des,
                         calmodel = ~region + utd + kjonn * ald_grp - 1)
tot_temp
tot_pop <- fill.template(universe = U, temp = tot_temp)
tot_pop
descal <- e.calibrate(design = des, df.population = tot_pop, calfun = "linear")
netto$kal_vekt <- weights(descal)
head(netto)
#
crossprod(netto$kal_vekt, netto$y1)
svystatTM(descal, y = ~y1, vartype = "cvpct")
sum(U$y1)
#
crossprod(netto$kal_vekt, netto$y2)
svystatTM(descal, y = ~y2, vartype = "cvpct")
sum(U$y2)
#
svystatTM(descal, y = ~y1, by = ~region, vartype = "cvpct")
tapply(U$y1, U$reg, sum)
#
svystatTM(descal, y = ~y, by = ~reg, estimator = "Mean", vartype = "cvpct")
tapply(d$y, d$reg, mean)
#
svystatTM(descal, y = ~y, by = ~edu, estimator = "Mean", vartype = "cvpct")
tapply(d$y, d$edu, mean)
#
svystatL(descal, expression(y / z), vartype = "cvpct")
sum(d$y) / sum(d$z)


  
