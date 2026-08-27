#-------------------------------------------------------
# R-kode for etterstratifisering hvis vi ikke har samme 
# trekksannsynlighet i alle etterstrata.                
# (Men det er enklere å bruke ReGenesees til dette også)
#-------------------------------------------------------
# Stratifisert tilfeldig trekking etter kjonn         
# krysset med aldersgruppene 16-24, 25-44, 45-66, 67+ 
# Total utvalgsstørrelse på 1000.                     
# Oversampling av de yngste og de eldste.              
#-------------------------------------------------------
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

HTestimator(utvalg$y2, utvalg$Prob) / 1e+6
sum(U$y2) / 1e+6

#------------------------
# Med frafall i utvalget 
#------------------------
# Antar at responssannsynligheten kun er avhengig av utdanningsnivå
utvalg$rs <- ifelse(utvalg$utd == "1", 0.30,
                    ifelse(utvalg$utd == "2", 0.40,
                           ifelse(utvalg$utd == "3", 0.50, 0.60)))
# Nettoutvalget
netto <- utvalg[runif(nrow(utvalg)) < utvalg$rs, ]
table(netto$kjonn, netto$ald_grp)

#-----------------------------------------
# Etterstratifisering: kjonn x utdanning 
#-----------------------------------------
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

#-----------------------------
# Etterstratifisering: region 
#-----------------------------
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
