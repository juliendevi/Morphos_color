# Chargement librairies 
library(pavo)
library(readxl)
library(factoextra)
library(rptR)

#### Jeu de données ####

# Importation des spectres de 300 à 700 nm
spectra = getspec("~/Documents/Stages/2021 (MNHN) - Coévolution vision-couleur des ailes chez les Morpho/spectro/aplomb/2021-03-30/", ext = "RFL8", lim = c(300,700))
Perou = subset(spectra, 'CL')
Guyane = subset(spectra, 'FG')


# Représenter les spectres en Overlay = tous sur mêmes axes x,y
plot(spectra, col=spec2rgb(spectra))

#### Smoothing spectra ####

#tests pour déterminer le span 
plotsmooth(spectra, minsmooth = 0.05, maxsmooth = 0.5, curves = 4, ask = FALSE)
#Smoothing
spectra_sm = procspec(spectra, opt = 'smooth', span = 0.2)
#plot
par(mfrow=c(1,1))
plot(spectra_sm, type='o', col=spec2rgb(spectra_sm))


##### Exploration par spécimen ####

# data brutes 
explorespec(spectra[,1:49], by=3, scale = "equal", legpos = "topright")
explorespec(spectra[,50:61], by=3, scale = "equal", legpos = "topright")


# data smooth
explorespec(spectra_sm, by=3, scale = "equal", legpos = "topright")


#### Répétabilité ####
infos <- read_excel("~/Documents/Stages/2021 (MNHN) - Coévolution vision-couleur des ailes chez les Morpho/spectro/aplomb/2021-03-30/2021-03-30.xlsx", sheet = "Feuil2")
test_repet = summary(spectra, subset = c('B1', 'S2', 'S3', 'H1'))
indv = infos$ID
esp = infos$sp
pop = infos$pop
test = cbind(test_repet, indv, esp)


is.factor(test$esp)
test$esp=as.factor(test$esp)

#### Boxplots
#Hue
ggplot(test, aes(x=esp, y=H1)) + geom_boxplot(width=0.5) + geom_jitter(aes(colour=indv), shape=16, size=2, width = 0) + xlab('Espèces') + ylab('Hue')
#Brightness
ggplot(test, aes(x=esp, y=B2)) + geom_boxplot(width=0.5) + geom_jitter(aes(colour=indv), shape=16, size=2, width = 0) + xlab('Espèces') + ylab('Brightness')
#Saturation
ggplot(test, aes(x=esp, y=S2)) + geom_boxplot(width=0.5) + geom_jitter(aes(colour=indv), shape=16, size=2, width = 0) + xlab('Espèces') + ylab('Saturation')
#Chroma
ggplot(test, aes(x=esp, y=S3)) + geom_boxplot(width=0.5) + geom_jitter(aes(colour=indv), shape=16, size=2, width = 0) + xlab('Espèces') + ylab('Chroma')

#### Répétabilité 
repH = rpt(H1~(1|indv), grname = 'indv', data=test, datatype = "Gaussian", nboot = 500)
repB = rpt(B1~(1|indv), grname = 'indv', data=test, datatype = "Gaussian", nboot = 500)
repSat = rpt(S2~(1|indv), grname = 'indv', data=test, datatype = "Gaussian", nboot = 500)
repC = rpt(S3~(1|indv), grname = 'indv', data=test, datatype = "Gaussian", nboot = 500)
repH
repB
repSat
repC


#### Calcul du coefficient de variation pour les 3 mesures ####

#H1
moyenne = by(test$H1, test$indv, mean)
ecart_type = by(test$H1, test$indv, sd)
coeffvar=ecart_type/moyenne
df = cbind(moyenne, ecart_type, coeffvar)
df
write.table(df,"~/Desktop/coeff_var_H1.txt", sep="\t")

#B1
moyenne = by(test$B1, test$indv, mean)
ecart_type = by(test$B1, test$indv, sd)
coeffvar=ecart_type/moyenne
df = cbind(moyenne, ecart_type, coeffvar)
df
write.table(df,"~/Desktop/coeff_var_B1.txt", sep="\t")

#S2
moyenne = by(test$S2, test$indv, mean)
ecart_type = by(test$S2, test$indv, sd)
coeffvar=ecart_type/moyenne
df = cbind(moyenne, ecart_type, coeffvar)
df
write.table(df,"~/Desktop/coeff_var_S2.txt", sep="\t")

#S3
moyenne = by(test$S3, test$indv, mean)
ecart_type = by(test$S3, test$indv, sd)
coeffvar=ecart_type/moyenne
df = cbind(moyenne, ecart_type, coeffvar)
df
write.table(df,"~/Desktop/coeff_var_S3.txt", sep="\t")

