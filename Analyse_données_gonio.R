#########################################################################################
########################### ANALYSES GONIO-SPECTROPHOTOMÈTRE ############################
#########################################################################################

##################### CHARGEMENT DES LIBRAIRIES #####################

library(tidyverse)
library(lightr)
library(pavo)
library(cowplot)

######### FONCTION RÉCUPÉRATION DES INFORMATIONS DU SPECTRE #########

# On va partitionner notre fichier .RAW8 pour récupérer les informations qui nous intéressent 
parse_RAW8 <- function(filename) {
  fichier = lr_parse_raw8(filename, specnum = 1) #Récupération spectre du 1er spectro (2ème == pondération fluctuations de la source)
  spectre = fichier[[1]] # [[2]] == métadonnées
  scope = spectre$scope # récupération du scope (nb photons)
  white = spectre$white # récupération du blanc 
  dark = spectre$dark # récupération du dark
  wl = spectre$wl
  data = data.frame(cbind(wl, scope, white, dark)) # tableau avec wl, scope, white & dark
  return(data) # retourne le tableau
}

########################### BRIGHTNESS ###########################

#### Récupérer les fichiers de mesure (hors blancs) ####

specfiles_B = list.files(path='.', pattern = "^[^H].[^WHITE].*\\.RAW8",full.names = TRUE, recursive = TRUE)
# pattern = ne commence pas par W et termine par .RAW8 files

#### Récupérer les angles de mesure ####

# Angle de la lumière incidente 
B_file_phi_inc = sapply(strsplit(specfiles_B, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-2]]))
# sapply = pour tous les éléments de la liste specfiles, et return result sous forme de vecteur
# strsplit = division de la chaine de caractère, dans specfiles, par la ponctuation ( '_' )
# Récupération sous forme numérique de l'avant dernier bout 

# Angle de la sonde collectrice 
B_file_phi_col = sapply(strsplit(specfiles_B, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-1]]))
# Récupération sous forme numérique du dernier bout

# calcul du tilt 
t = (B_file_phi_inc[4] - B_file_phi_col[4])/2

# Calcul de la normale lors des prises de mesure 
B_file_norm = (B_file_phi_inc - B_file_phi_col)/2 - t 

#### Récupérer d'autres infos dans le titre du fichier ####

# Ajouter pour sortir ID, nom espèce et sexe ? 
#file_ID = sapply(strsplit(specfiles_B, "[[:punct:]]"),function(x) as.character(x[[length(x)-3]]))


preprocess_renorm_B <- function(file) {
  file_infos = strsplit(file, "[[:punct:]]")[[1]] # On split les infos contenues dans le nom de fichier
  norm = (as.numeric(file_infos[length(file_infos)-2]) - as.numeric(file_infos[length(file_infos)-1]))/2 - t # On calcul la normale à partir des angles d'incidence et de collecte
  whitefiles = list.files(path = ".",pattern = "^B_WHITE", full.names = TRUE) # On récupère la liste des blancs
  white_phi_inc = sapply(strsplit(whitefiles, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-2]])) # Récupération de leur angle d'incidence
  white_phi_col = sapply(strsplit(whitefiles, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-1]])) # Récupération de leur angle de collecte
  white_norms = (white_phi_inc - white_phi_col)/2 # Calcul de la normale pour chaque blanc 
  spectre_df = parse_RAW8(file) # On récupère et partitionne le spectre à corriger
  white_df = parse_RAW8(whitefiles[which(white_norms==norm)]) # Choix du bon blanc à appliquer 
  correction = (spectre_df$scope - spectre_df$dark)/(white_df$scope - spectre_df$dark)*100 # Application de la correction au spectre
  spectre_final = as.rspec(cbind(spectre_df$wl, correction), lim=c(300,700)) # Enregistrer sous forme de spectre avec wl & processed 
  write.csv(spectre_final, gsub("\\.RAW8$", ".csv", file),row.names = FALSE) # écrire un csv
}

#### Application de la correction à tous les spectres de brightness ####

sapply(specfiles_B, preprocess_renorm_B) # Application de la correction à chaque item de la liste

########################### HUE ###########################

#### Récupérer les fichiers de mesure (hors blancs) ####

specfiles_H = list.files(path = ".", pattern = "^[^B].[^WHITE].*\\.RAW8",full.names = TRUE, recursive = TRUE)

#### Récupérer les angles de mesure ####

# Angle de la lumière incidente 
H_file_phi_inc = sapply(strsplit(specfiles_H, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-2]]))
# sapply = pour tous les éléments de la liste specfiles, et return result sous forme de vecteur
# strsplit = division de la chaine de caractère, dans specfiles, par la ponctuation ( '_' )
# Récupération sous forme numérique de l'avant dernier bout 

# Angle de la sonde collectrice 
H_file_phi_col = sapply(strsplit(specfiles_H, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-1]]))
# Récupération sous forme numérique du dernier bout

# Calcul de la normale lors des prises de mesure 
H_file_span = (H_file_phi_col + H_file_phi_inc) + 180 - 720

#### Fonction de renormalisation des spectres de Brightness ####

preprocess_renorm_H <- function(file) {
  file_infos = strsplit(file, "[[:punct:]]")[[1]] # On split les infos contenues dans le nom de fichier
  file_span = as.numeric(file_infos[length(file_infos)-2]) + as.numeric(file_infos[length(file_infos)-1]) + 180 - 720
  whitefiles = list.files(path = ".", pattern = "^H_WHITE", full.names = TRUE) # Récupération des blancs
  white_phi_inc = sapply(strsplit(whitefiles, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-2]])) # Récupération de leur angle d'incidence
  white_phi_col = sapply(strsplit(whitefiles, "[[:punct:]]"),function(x) as.numeric(x[[length(x)-1]])) # Récupération de leur angle de collecte
  white_entre_axe = white_phi_col + white_phi_inc + 180 - 720 # Calcul de l'entre-axe
  spectre_df = parse_RAW8(file) # On récupère et partitionne le spectre à corriger
  white_df = parse_RAW8(whitefiles[which(white_entre_axe==file_span)]) # Coix du bon blanc pour correction
  correction = (spectre_df$scope - spectre_df$dark)/(white_df$scope - spectre_df$dark)*100 # Application de la correction au spectre
  spectre_final = as.rspec(cbind(spectre_df$wl, correction), lim=c(300,700)) # Enregistrer sous forme de spectre avec wl & processed 
  write.csv(spectre_final, gsub("\\.RAW8$", ".csv", file),row.names = FALSE) # écrire un csv
}

#### Application de la correction à tous les spectres de brightness ####

sapply(specfiles_H, preprocess_renorm_H) # Application de la correction à chque item de la liste

##################### ESTIMATION DES PARAMÈTRES D'IRIDESCENCE #####################

#### Extraire les variables H1 et B2 ####

get_col_var <- function(variable){
  spectre = suppressWarnings(getspec(sep = ",", ext = "csv", subdir = TRUE, subdir.names = FALSE)) # Récupération de tous les spectres 
  spectre_variable = subset(spectre, variable) # subset variable qui nous intéresse (H or B)
  spectre_variable = procspec(spectre_variable, 'smooth', 'zero') # lissage + min=0
  col_var = summary(spectre_variable, subset = c("H1", "B2"), wlmin=350, wlmax=700) # Récupération des variables 
  col_var$I = sapply(strsplit(rownames(col_var), "_"),function(x) as.numeric(x[length(x)-1]))
  col_var$C = sapply(strsplit(rownames(col_var), "_"),function(x) as.numeric(x[length(x)]))
  col_var$span = 180 + col_var$I + col_var$C - 720
  col_var$halfspan = col_var$span / 2
  col_var$normale = (col_var$I - col_var$C) / 2
  # Si B2 trop faible = noir donc paramètre B non applicable 
  col_var[col_var$B2<8.5, "H1"] = NA
  # Remove artefact at the edges
  col_var$H1[col_var$H1 %in% c(300,700)] = NA
  # Discard large angles
  col_var = col_var[col_var$span<90,]
  return(col_var)
}

#### Fonctions normale + cosinus pour l'estimation des paramètres ####

fnorm = function(x, Bmax, t, gammaB) {
  Bmax * exp(-0.5*(x-t)^2/gammaB^2)
}

fcos = function(x, Hmax, gammaH) {
  Hmax * cos(gammaH * x / 180 * pi)
}

#### Fonction d'estimation optimale des paramètres pour B et H par méthode non linéaire des moindres carrés ####

# Pour Brightness
find_params_nls_normale = function(brightness_folder) {
  maxi = max(brightness_folder$B2)
  norm = brightness_folder$normale[which.max(brightness_folder$B2)]
  sigm = abs(brightness_folder$normale[which.min(abs(brightness_folder$B2 - exp(-0.5) * maxi))] - norm)
  fit = nls(B2 ~ fnorm(normale, Bmax, t, gammaB),
            data = brightness_folder,
            start = c("Bmax"=maxi, "t"=norm, "gammaB"=sigm),
            lower = c("Bmax"=0 , "t"=-50 , "gammaB"=0),
            algorithm = "port",
            nls.control(warnOnly = TRUE))
  return(summary(fit)$coefficients[,1])
}

# Pour Hue 
find_params_nls_span = function(hue_folder) {
  hue_folder = hue_folder[!is.na(hue_folder$H1),]
  if (nrow(hue_folder)<2) {
    # If only one measurement, we can't estimate parameters
    return(rep(NA,2))
  } else {
    maxi = max(hue_folder$H1)
    s = 0.6
    fit = nls(H1 ~ fcos(halfspan, Hmax, gammaH),
              data = hue_folder,
              start = c("Hmax"=maxi, "gammaH"=s),
              control = nls.control(warnOnly = TRUE))
    return(summary(fit)$coefficients[,1])
  }
}


##################### FONCTION PLOT VIA GGPLOT2 #####################

ggplot_rspec = function(x) {
  rspecdata = as.rspec(x)
  rspecdata_long = tidyr::gather(rspecdata, name, spec, -wl)
  g = ggplot(data = rspecdata_long,
             aes(x = wl, y = spec, group = name, col = name)) +
    geom_line() +
    ylab("Reflectance (%)") +
    xlab("Wavelength") +
    scale_colour_manual(values = spec2rgb(rspecdata))
  return(g)
}


##################### VISUALISATION DES SPECTRES #####################

theme_set(theme_minimal())

#### Récupération des spectres corrigés ####
spec = getspec(sep=',', ext = 'csv')
hue_spec = subset(spec, 'H')
brightness_spec = subset(spec, 'B')

#### Lissage des spectres et valeur min = 0 ####
sm_hue_spec = procspec(hue_spec, 'smooth', 'zero')
sm_brightness_spec = procspec(brightness_spec, 'smooth', 'zero')

#### Plot ####

ggplot_rspec(sm_hue_spec) + theme(legend.position = "none")
ggplot_rspec(sm_brightness_spec) + theme(legend.position = "none")

#### On récupère les variables ####

hue_var = get_col_var('H') # variables Hue
hue_var$rgb = spec2rgb(sm_hue_spec) # application code couleur

brightness_var = get_col_var('B')  # variables Brightness
brightness_var$rgb = spec2rgb(sm_brightness_spec) # application code couleur


#### On réalise les régressions ####

Hue_reg = find_params_nls_span(hue_var)
Brightness_reg = find_params_nls_normale(brightness_var)

#### Plot for Hue ####

# Estimation variation de la B selon les angles
vB = ggplot(hue_var, aes(x = halfspan, y = B2, col = factor(halfspan))) +
  geom_point(size = 3) +
  xlab(expression((Phi[inc]-Phi[col])/2)) +
  ylim(c(0, max(brightness_var$B2, hue_var$B2))) +
  ylab("Brightness B (%)") +
  #scale_color_manual(values = unname(spec2rgb(sm_hue_spec))) +
  theme(legend.position = "none")

# Variation Hue 
H=ggplot(hue_var, aes(x = halfspan, y = H1, col = factor(halfspan))) +
  xlab(expression((Phi[inc]-Phi[col])/2)) +
  stat_function(fun = fcos, args = Hue_reg, color = "red") +
  geom_point(size = 3) +
  ylim(range(brightness_var$H1, hue_var$H1, na.rm = TRUE)) +
  ylab("Hue H (nm)") +
  scale_color_manual(values = unname(spec2rgb(sm_hue_spec))) +
  theme(legend.position = "none") +
  annotate("text", x = 11, y = 460, color = "red",
           label = sprintf("H[max]== %.0f~nm", Hue_reg[["Hmax"]]),
           parse = TRUE) +
  annotate("text", x = 11, y = 455, color = "red",
           label = sprintf("gamma[H]== %.2f", Hue_reg[["gammaH"]]),
           parse = TRUE)

 #### Plot for Brightness ####

B = ggplot(brightness_var, aes(x = normale, y = B2, col = factor(normale))) +
  stat_function(fun = fnorm, args = Brightness_reg, color = "red") +
  geom_point(shape = "square", size = 3) +
  xlab(expression((Phi[inc]+Phi[col])/2)) +
  ylim(c(0, max(brightness_var$B2, hue_var$B2))) +
  ylab("Brightness B (%)") +
  scale_color_manual(values = unname(spec2rgb(sm_brightness_spec))) +
  theme(legend.position = "none") +
  annotate("text", x = 5, y = 80, color = "red",
           label = sprintf("B[max]== %.0f*'%%'", Brightness_reg[["Bmax"]]),
           parse = TRUE) +
  annotate("text", x = 5, y = 70, color = "red",
           label = sprintf("gamma[B]== %.2f", Brightness_reg[["gammaB"]]),
           parse = TRUE) +
  annotate("text", x = 5, y = 60, color = "red",
           label = sprintf("t== %.0f*'°'", Brightness_reg[["t"]]), parse = TRUE)

# Estimation variation Hue avec l'angle 
vH = ggplot(brightness_var, aes(x = normale, y = H1, col = factor(normale))) +
  geom_point(shape = "square", size = 3) +
  xlab(expression((Phi[inc]+Phi[col])/2)) +
  ylim(range(brightness_var$H1, hue_var$H1, na.rm = TRUE)) +
  ylab("Hue H (nm)") +
  #scale_color_manual(values = unname(spec2rgb(sm_brightness_spec))) +
  theme(legend.position = "none")

# Plot des 4 graphes d'estimations des paramètres 
plot_grid(H, vB, B, vH, labels=c('A', 'B', 'C', 'D'), ncol=2, nrow = 2)



