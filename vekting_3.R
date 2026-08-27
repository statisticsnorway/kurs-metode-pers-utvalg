#--------------------------------------------------------------------------
# Kalibrering med ReGenesees. Kalibreringsmodell: region + kjonn * ald_grp 
#--------------------------------------------------------------------------
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

summary(netto$d_vekt / mean(netto$d_vekt) * mean(netto$kal_vekt) / netto$kal_vekt)

# Estimat for hele populasjonen
svystatTM(descal, y = ~y1, estimator = "Mean", vartype = "se")
mean(U$y1)
# sjekk også weighted.mean(netto$y1, w = netto$kal_vekt)

# Estimat etter region
svystatTM(descal, y = ~y1, by = ~region, estimator = "Mean", vartype = "se")
tapply(U$y1, U$region, mean)

# Estimat etter utdanning
svystatTM(descal, y = ~y1, by = ~utd, estimator = "Mean", vartype = "se")
tapply(U$y1, U$utd, mean)

# Estimat etter kjønn krysset med aldersgruppe
svystatTM(descal, y = ~y1, by = ~kjonn * ald_grp, estimator = "Mean", vartype = "se")
tapply(U$y1, list(U$kjonn, U$ald_grp), mean)

# Estimat - forholdet mellom to populasjonssummer/-gjennomsnitt
svystatL(descal, expression(y2 / y3), vartype = "se")
sum(U$y2) / sum(U$y3)

svystatL(descal, expression(y2 / y3), by = ~region, vartype = "se", conf.int = TRUE)
tapply(U$y2, U$region, sum) / tapply(U$y3, U$region, sum)

#------------------------------------------------------------
# Kalibrering med ReGenesees: region + utd + kjonn * ald_grp 
#------------------------------------------------------------
des <- e.svydesign(data = netto, ids = ~pnr, weight = ~d_vekt)
tot_temp <- pop.template(data = des,
                         calmodel = ~region + utd + kjonn * ald_grp - 1)
tot_temp
tot_pop <- fill.template(universe = U, temp = tot_temp)
tot_pop
descal <- e.calibrate(design = des, df.population = tot_pop, calfun = "linear")
netto$kal2_vekt <- weights(descal)
head(netto)

summary(netto$d_vekt / mean(netto$d_vekt) * mean(netto$kal2_vekt) / netto$kal2_vekt)

svystatTM(descal, y = ~y1, estimator = "Mean", vartype = "se")
mean(U$y1)

svystatTM(descal, y = ~y1, by = ~region, estimator = "Mean", vartype = "se")
tapply(U$y1, U$region, mean)

svystatTM(descal, y = ~y1, by = ~utd, estimator = "Mean", vartype = "se")
tapply(U$y1, U$utd, mean)

svystatL(descal, expression(y2 / y3), vartype = "se")
sum(U$y2) / sum(U$y3)

#--------------------------------------------
# Dobbeltkalibrering med ReGenesees:                       
# 1) Frafallsvekting mot utdanning.         
# 2) Kalibrering mot region + kjonn x alder.
#--------------------------------------------
des <- e.svydesign(data = netto, ids = ~pnr, weight = ~d_vekt)
tot_temp <- pop.template(data = des,
                         calmodel = ~utd - 1)
tot_temp
tot_utv <- fill.template(universe = utvalg, temp = tot_temp)
tot_utv
descal_1 <- e.calibrate(design = des, df.population = tot_utv, calfun = "linear")

tot_temp <- pop.template(data = descal_1,
                         calmodel = ~region + kjonn * ald_grp - 1)
tot_temp
tot_pop <- fill.template(universe = U, temp = tot_temp)
tot_pop
descal_2 <- e.calibrate(design = descal_1, df.population = tot_pop, calfun = "linear")
netto$kal3_vekt <- weights(descal_2)
head(netto)

summary(netto$d_vekt / mean(netto$d_vekt) * mean(netto$kal3_vekt) / netto$kal3_vekt)

svystatTM(descal_2, y = ~y1, estimator = "Mean", vartype = "se")
mean(U$y1)

svystatTM(descal_2, y = ~y1, by = ~region, estimator = "Mean", vartype = "se")
tapply(U$y1, U$region, mean)

svystatTM(descal_2, y = ~y1, by = ~utd, estimator = "Mean", vartype = "se")
tapply(U$y1, U$utd, mean)

svystatL(descal, expression(y2 / y3), vartype = "se")
sum(U$y2) / sum(U$y3)








