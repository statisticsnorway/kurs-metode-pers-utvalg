#--------------------------------
# Etterstratifisering mot region 
#--------------------------------
t_U <- table(U$region)
t_n <- table(netto$region)
w <- t_U / t_n
w
netto$es_vekt <- w[netto$region]

# Estimering av gjennomsnitt
weighted.mean(netto$y1, w = netto$es_vekt)
mean(U$y1)

weighted.mean(netto$y2, w = netto$es_vekt)
mean(U$y2)

weighted.mean(netto$y3, w = netto$es_vekt)
mean(U$y3)

#----------------------------------------------
# Etterstratifisering mot kjonn x aldersgruppe 
#----------------------------------------------
t_U <- table(U$kjonn, U$ald_grp)
t_n <- table(netto$kjonn, netto$ald_grp)
w <- t_U / t_n
w
# Den etterstratifiserte vekten
netto$es_vekt <- w[cbind(netto$kjonn, netto$ald_grp)]

# Estimering av gjennomsnitt
weighted.mean(netto$y1, w = netto$es_vekt)
mean(U$y1)

weighted.mean(netto$y2, w = netto$es_vekt)
mean(U$y2)

weighted.mean(netto$y3, w = netto$es_vekt)
mean(U$y3)

#---------------------------------------
# Etterstratifisering: kjonn x utdanning 
#---------------------------------------
t_U <- table(U$kjonn, U$utd)
t_n <- table(netto$kjonn, netto$utd)
w <- t_U / t_n
w
netto$es_vekt <- w[cbind(netto$kjonn, netto$utd)]

# Estimering av gjennomsnitt
weighted.mean(netto$y1, w = netto$es_vekt)
mean(U$y1)

weighted.mean(netto$y2, w = netto$es_vekt)
mean(U$y2)

weighted.mean(netto$y3, w = netto$es_vekt)
mean(U$y3)

