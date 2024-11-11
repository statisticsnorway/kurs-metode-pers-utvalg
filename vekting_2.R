
##########################################
# Etterstratifisering: kjonn x utdanning #
##########################################
t_U <- table(U$kjonn, U$utd)
t_n <- table(netto$kjonn, netto$utd)
w <- t_U / t_n
w
netto$es_vekt <- w[cbind(netto$kjonn, netto$utd)]
#
# Estimering av populasjonstotaler (vektet sum)
crossprod(netto$es_vekt, netto$y1) / 1e+6
sum(U$y1) / 1e+6
#
crossprod(netto$es_vekt, netto$y2) / 1e+6
sum(U$y2) / 1e+6

##############################################################
# Kalibrering med ReGenesees: region + utd + kjonn * ald_grp #
##############################################################
des <- e.svydesign(data = netto, ids = ~pnr, weight = ~d_vekt)
tot_temp <- pop.template(data = des,
                         calmodel = ~region + utd + kjonn * ald_grp - 1)
tot_temp
tot_pop <- fill.template(universe = U, temp = tot_temp)
tot_pop
descal <- e.calibrate(design = des, df.population = tot_pop, calfun = "linear")
netto$kal_vekt <- weights(descal)
head(netto)
#
summary(netto$d_vekt / mean(netto$d_vekt) * mean(netto$kal_vekt) / netto$kal_vekt)
#
#crossprod(netto$kal_vekt, netto$y1)
svystatTM(descal, y = ~y1, vartype = "se") / 1e+6
sum(U$y1) / 1e+6
#
#crossprod(netto$kal_vekt, netto$y2)
svystatTM(descal, y = ~y2, vartype = "se") / 1e+6
sum(U$y2) / 1e+6
#
svystatTM(descal, y = ~y1, by = ~region, vartype = "se")
tapply(U$y1, U$reg, sum)
#
svystatTM(descal, y = ~y1, by = ~region, estimator = "Mean", vartype = "se")
tapply(U$y1, U$region, mean)
#
svystatTM(descal, y = ~y1, by = ~utd, estimator = "Mean", vartype = "se")
tapply(U$y1, U$utd, mean)
#
svystatL(descal, expression(y1 / y2), vartype = "se")
sum(U$y1) / sum(U$z)

##############################################
# Dobbelkalibrering:                         #
# 1) Frafallsvekting mot utdanning.          #
# 2) Kalibrering mot region + kjonn x alder. #
##############################################
des <- e.svydesign(data = netto, ids = ~pnr, weight = ~d_vekt)
tot_temp <- pop.template(data = des,
                         calmodel = ~utd - 1)
tot_temp
tot_utv <- fill.template(universe = utvalg, temp = tot_temp)
tot_utv
descal_1 <- e.calibrate(design = des, df.population = tot_utv, calfun = "linear")
#
tot_temp <- pop.template(data = descal_1,
                         calmodel = ~region + kjonn * ald_grp - 1)
tot_temp
tot_pop <- fill.template(universe = U, temp = tot_temp)
tot_pop
descal_2 <- e.calibrate(design = descal_1, df.population = tot_pop, calfun = "linear")
netto$kal2_vekt <- weights(descal_2)
head(netto)
#
summary(netto$d_vekt / mean(netto$d_vekt) * mean(netto$kal2_vekt) / netto$kal2_vekt)
#
svystatTM(descal_2, y = ~y1, vartype = "se") / 1e+6
sum(U$y1) / 1e+6
#
svystatTM(descal, y = ~y2, vartype = "se") / 1e+6
sum(U$y2) / 1e+6
#
svystatTM(descal, y = ~y1, by = ~region, vartype = "se")
tapply(U$y1, U$reg, sum)
#
svystatTM(descal, y = ~y1, by = ~region, estimator = "Mean", vartype = "se")
tapply(U$y1, U$region, mean)
#
svystatTM(descal, y = ~y1, by = ~utd, estimator = "Mean", vartype = "se")
tapply(U$y1, U$utd, mean)
#
svystatL(descal, expression(y1 / y2), vartype = "se")
sum(U$y1) / sum(U$z)


#######################################################
# R-kode for etterstratifisering hvis vi ikke har lik #
# trekksannsynlighet innen hvert etterstrata          #
#######################################################
# Stratifisert tilfeldig trekking etter kjonn         #
# krysset med aldersgruppene 16-24, 25-44, 45-66, 67+ #
# Total utvalgsstørrelse på 1000.                     #
# Oversampling av de yngste og de eldste.             # 
#######################################################
U <- d[d$alder >= 16, ]
U$ald_grp <- cut(U$alder, breaks = c(15, 24, 44, 66, 100))
table(U$kjonn, U$ald_grp)
n <- c(150, 100, 100, 150,
       150, 100, 100, 150)
U <- U[order(U$kjonn, U$ald_grp), ]
s <- strata(data = U,
            stratanames = c("kjonn", "ald_grp"),
            size = n,
            method = "srswor")
# Legger på alle variabler fra trekkerammen
utvalg <- getdata(U, s)
# Designvekten
utvalg$d_vekt <- 1 / utvalg$Prob
head(utvalg)
# Estimerer populasjonstotal (i millioner) med Horvitz-Thomsen-estimatoren
HTestimator(utvalg$y1, utvalg$Prob) / 1e+6
# Den faktiske totalen er:
sum(U$y1) / 1e+6
#
HTestimator(utvalg$y2, utvalg$Prob) / 1e+6
sum(U$y2) / 1e+6

##########################
# Med frafall i utvalget #
##########################
# Antar at responssannsynligheten kun er avhengig av utdanningsnivå
utvalg$rs <- ifelse(utvalg$utd == "1", 0.30,
                    ifelse(utvalg$utd == "2", 0.40,
                           ifelse(utvalg$utd == "3", 0.50, 0.60)))
# Nettoutvalget
netto <- utvalg[runif(nrow(utvalg)) < utvalg$rs, ]
table(netto$kjonn, netto$ald_grp)

##########################################
# Etterstratifisering: kjonn x utdanning #
##########################################
t_U <- table(U$kjonn, U$utd)
t_n <- tapply(netto$d_vekt, list(netto$kjonn, netto$utd), sum)
w <- t_U / t_n
w
# Etterstratifisert vekt
netto$es_vekt <- netto$d_vekt * w[cbind(netto$kjonn, netto$utd)]
#
# Estimering av populasjonstotaler (vektet sum)
crossprod(netto$es_vekt, netto$y1) / 1e+6
sum(U$y1) / 1e+6
#
crossprod(netto$es_vekt, netto$y2) / 1e+6
sum(U$y2) / 1e+6

###############################
# Etterstratifisering: region #
###############################
t_U <- table(U$region)
t_n <- tapply(netto$d_vekt, netto$region, sum)
w <- t_U / t_n
w
# Etterstratifisert vekt
netto$es_vekt <- netto$d_vekt * w[netto$region]
#
# Estimering av populasjonstotaler (vektet sum)
crossprod(netto$es_vekt, netto$y1) / 1e+6
sum(U$y1) / 1e+6
#
crossprod(netto$es_vekt, netto$y2) / 1e+6
sum(U$y2) / 1e+6




