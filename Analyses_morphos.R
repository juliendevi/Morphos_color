##############################################################################
###################### ANALYSES DES SPECTRES À L'APLOMB ######################
##############################################################################

######################## Chargement librairies ###############################
library(pavo)
library(factoextra)
library(tidyverse)
library(phytools)
library(picante)
library(geiger)
library(geomorph)
library(ggtree)
library(cowplot)
library(RColorBrewer)
library(ggphylomorpho)
library(lemon)
library(ggimage)

######################## Chargement fonctions ###############################

# Test de Kmult d'Adams 
source("~/Desktop/tuto/Kmult adams.R") 
# Fonction de calcul de moyennes multivariées
source("~/Desktop/tuto/meangr.R") # Vincent Debat #

############### Theme minimal + set-working directory ########################
theme_set(theme_minimal())
setwd('~/Desktop/Macroevol/Cu1_M3')

################## Fonction pour afficher les spectres #####################

## Modifié d'après Gruson et al., 2018 ## Hugo Gruson ##
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

###### Fonction pour récupérer les informations de chaque individu ########

recup_infos = function(specfile) {
  file_infos = strsplit(specfile, "[/_\\]")[[1]]
  canop = list(c('theseus', 'niepelti', 'amphitryon', 'telemachus', 'hercules', 'cisseis', 'hecuba', 'anaxibia', 'cypris', 'rhetenor'))
  ID = file_infos[2]
  espece = file_infos[3]
  sex = file_infos[4]
  couleur = file_infos[5]
  if (espece %in% canop[[1]]){
    habitat='Canopy'
  }else{
    habitat = 'Understory'
  }
  data.frame(ID, espece, sex, habitat, couleur)
}

###### Transposition des spectres pour analyse et fonction inverse ######

spec_to_table <- function(x){
  table = t(x)
  colnames(table)=table[1,]
  table = table[-1,]
  return(table)
}

table_to_spec <- function(x){
  spec = t(x)
  wl = as.numeric(rownames(spec))
  spec = cbind(wl, spec)
  rownames(spec)=NULL
  return(spec)
}

########################## Arbre phylogénétique ###########################

tree_morphos = read.nexus("~/Desktop/tuto/morpho_tree2014bis.nex") # Chazot et al., 2016 #

##################### Jeu de données macroévolution ######################

# Localisation des fichiers 
specfiles = list.files(path='.', pattern='*\\.RFL8',full.names = TRUE, recursive = TRUE)
# Récupération des infos 
infos = do.call(rbind, lapply(specfiles, recup_infos))
# Importation des spectres de 300 à 700 nm
spectra = getspec("~/Desktop/Macroevol/Cu1_M3", ext = "RFL8", lim = c(300,700), subdir = T)
# Smooth de tous les spectres avec le même coefficient 
spectra_sm = procspec(spectra, opt = 'smooth', fixneg = 'zero', span = 0.2)

# Représentation graphique de l'ensemble des 277 spectres (278 si supplémentaire ajouté)
plot(spectra_sm)

#### Base de Figure box1 ####
a=c(1,5)
c = summary(spectra_sm, subset=c('H1', 'B1', 'S2', 'S3'))
c = c[4,]
Rmin = min(spectra_sm$BC084_aega_Male_bleu)
ggplot_rspec(spectra_sm[,a])+
  theme_linedraw()+theme(legend.position = 'none') + 
  geom_area(fill="lightgrey") + 
  geom_vline(xintercept = c$H1, linetype="dashed") + 
  geom_vline(xintercept = c$H1-50, linetype="dashed") + 
  geom_vline(xintercept = c$H1+50, linetype="dashed") + 
  geom_hline(yintercept = spectra_sm[which(spectra_sm$wl==c$H1),5], linetype="dashed") + 
  geom_hline(yintercept = spectra_sm[which(spectra_sm$BC084_aega_Male_bleu==Rmin),5], linetype="dashed")

################## MANNOVA non phylogénétique ALL #########################

prepACP = procspec(spectra_sm, opt = c('bin', 'center')) 
# Transposition wavelength en variables
prepACP = spec_to_table(prepACP)

# Matrice de corrélation 
cor(prepACP)

# MANOVA effet sex + habitat (C vs U)
summary(manova(prepACP~infos$sex))
summary(manova(prepACP~infos$habitat))
summary(manova(prepACP~infos$sex*infos$habitat))
summary(manova(prepACP~infos$espece))
summary(manova(prepACP~infos$espece*infos$sex))
# Différence significative entre sexes et espèces => besoin de séparer les 2 sexes

############### Calcul spectres moyennes par espèce ALL ######################

# Calcul des spectres moyens 
spectres_mean = aggspec(spectra_sm, by=infos$espece, FUN = mean)
# Représentation graphique 
ggplot_rspec(spectres_mean) #+ ggtitle('spectres moyens de chaque espèces')

#################### Sous jeu de données Males #######################

# Sous échantillonnage des spectres 
spectres_males = subset(spectra_sm, 'Male')
# Récupération des infos 
specfiles_Male = list.files(path='.', pattern='Male.*\\.RFL8',full.names = TRUE, recursive = TRUE)
# Récupération des infos 
infos_males = do.call(rbind, lapply(specfiles_Male, recup_infos))

############### Calcul spectres moyennes par espèce (males) ######################

# Calcul des spectres moyens 
spectres_mean_males = aggspec(spectres_males, by=infos_males$espece, FUN = mean)
# Représentation graphique 
library(directlabels)
a = ggplot_rspec(spectres_mean_males) #+ ggtitle('Mean spectres des males de chaque espèce')
direct.label(a)

######## Reconstruction de l'arbre phylogénétique avec la couleur MALE ###############

## ATTENTION ! la plupart des fonctions nécessitent variables en colonnes et individus en lignes
# C'est l'inverse des spectres (cf. modifications pour ACP)

# Transposition du dataframe
mean_males = spec_to_table(spectres_mean_males)

## ATTENTION ! il faut vérifier que les labels de l'arbre et des data sont dans le même sens !

# Mise dans l'ordre des spectres moyens
couleur_males = match.phylo.data(tree_morphos, mean_males)$data
# Mise dans l'ordre des labels de l'arbre
tree = match.phylo.data(tree_morphos, mean_males)$phy


# Spectres moyens remis sous forme rspec
couleur_males = table_to_spec(couleur_males)

# Visualisation de l'arbre et des couelurs dans Figure 1 (fin de section)

#################### ACP spectre moyen par espèce Males #########################

# ACP
ACP_mean_males = prcomp(mean_males, scale. = T)
# Résumé de l'ACP
summary(ACP_mean_males)

### Graphes ACP 
# Explication des axes 
fviz_screeplot(ACP_mean_males, addlabels = T) + ggtitle('')
# Représentation individus
fviz_pca_ind(ACP_mean_males, asp=1, pointshape=20) + scale_color_discrete() + ggtitle('ACP mean males')
fviz_pca_ind(ACP_mean_males, axes = c(2,3), asp=1, pointshape=20) + ggtitle('ACP mean males')

######################### Phylomorphospace Males #############################

# Mise en couleur selon type de vol 
mypalette<-c("#357AB7", "#82C46C")
# Vérification des noms données vs. arbre phylogénétique 
name.check(tree_morphos, mean_males)
phylo = as.data.frame(cbind(ACP_mean_males$x[,1:2]))

for (i in 1:nrow(phylo)){
  canop = list(c('theseus', 'niepelti', 'amphitryon', 'telemachus', 'hercules', 'cisseis', 'hecuba', 'anaxibia', 'cypris', 'rhetenor'))
  if (row.names(phylo[i,]) %in% canop[[1]] ){
    phylo$habitat[i]='Canopée'
  }
  else {
    phylo$habitat[i]='Sous-bois'
  }
}

phylomorphomales = ggphylomorpho(tree_morphos, phylo, xvar = phylo$PC1, yvar=phylo$PC2, factorvar = phylo$habitat, labelvar = row.names(phylo), edge.width = 0.5, title = 'a) Mâles', xlab = 'PC1 (61,2%)', ylab = 'PC2 (24,7%)') + scale_color_manual(values=mypalette, name='Habitat') + theme(legend.position = 'none', text = element_text(size=14))

######################### Signal phylogénétique Males #############################

Test.Kmult(mean_males, tree_morphos, iter=999)
# Y a du signal phylogénétique (1,26 = fort)
# Besoin prendre en compte la phylogénie dans nos analyses => MANNOVA phylogénétique 

######################### MANNOVA phylogénétique Males #############################
a=as.matrix(ACP_mean_males$x[,1:7])
voltypr=as.factor(phylo$habitat)
names(voltypr)<-dimnames(a)[[1]]
aov.phylo(a~voltypr, tree_morphos, nsim = 10, test="Wilks") 
# MANOVA non phylo
summary(manova(ACP_mean_males$x[,1:7]~voltypr))

#################### Sous jeu de données Femelles #######################

# Sous échantillonnage des spectres 
spectres_femelles = subset(spectra_sm, 'Femelle')
# Récupération des infos 
specfiles_Femelles = list.files(path='.', pattern='Femelle.*\\.RFL8',full.names = TRUE, recursive = TRUE)
# Récupération des infos 
infos_femelles = do.call(rbind, lapply(specfiles_Femelles, recup_infos))

############### Calcul spectres moyennes par espèce (femelles) ######################

# Calcul des spectres moyens 
spectres_mean_femelles = aggspec(spectres_femelles, by=infos_femelles$espece, FUN = mean)
# Représentation graphique 
ggplot_rspec(spectres_mean_femelles) + ggtitle('Mean spectres femelles par espèce')

######## Reconstruction de l'arbre phylogénétique avec la couleur FEMELLE ###############

## ATTENTION ! la plupart des fonctions nécessitent variables en colonnes et individus en lignes
# C'est l'inverse des spectres (cf. modifications pour ACP)

# Transposition du dataframe
mean_femelles = spec_to_table(spectres_mean_femelles)

## ATTENTION ! il faut vérifier que les labels de l'arbre et des data sont dans le même sens !

# Mise dans l'ordre des spectres moyens
couleur_femelles = match.phylo.data(tree_morphos, mean_femelles)$data
# Mise dans l'ordre des labels de l'arbre
tree = match.phylo.data(tree_morphos, mean_femelles)$phy

# Spectres moyens remis sous forme rspec
couleur_femelles = table_to_spec(couleur_femelles)

# On affiche l'arbre 
plot(tree,"p", FALSE, font = 1, label.offset = 3, x.lim = 35, no.margin = FALSE)
# On ajoute les estimations de couleur 
tiplabels(pch=22, cex=3, bg=spec2rgb(as.rspec(couleur_femelles)),  adj = 1.4)
title(main="Phylogénie approximation des spectres couleur dorsale des Femelles", font.main=3)

#################### ACP spectre moyen par espèce Femelles #########################

# ACP
ACP_mean_femelles = prcomp(mean_femelles, scale. = T)
# Résumé de l'ACP
summary(ACP_mean_femelles)

### Graphes ACP 
# Explication des axes 
fviz_screeplot(ACP_mean_femelles, addlabels = T) 
# Représentation individus
fviz_pca_ind(ACP_mean_femelles, asp=1, pointshape=20) + scale_color_discrete() + ggtitle('ACP mean femelles')
fviz_pca_ind(ACP_mean_femelles, axes = c(2,3), asp=1, pointshape=20) + ggtitle('ACP mean femelles')


######################### Phylomorphospece Femelles #############################

# Vérification des noms données vs. arbre phylogénétique 
name.check(tree_morphos, mean_femelles)
phylo_fem = as.data.frame(cbind(ACP_mean_femelles$x[,1:2]))

for (i in 1:nrow(phylo_fem)){
  canop = list(c('theseus', 'niepelti', 'amphitryon', 'telemachus', 'hercules', 'cisseis', 'hecuba', 'anaxibia', 'cypris', 'rhetenor'))
  if (row.names(phylo_fem[i,]) %in% canop[[1]] ){
    phylo_fem$habitat[i]='Canopée'
  }
  else {
    phylo_fem$habitat[i]='Sous-bois'
  }
}

phylomorphofemelle = ggphylomorpho(tree_morphos, phylo_fem, xvar = phylo_fem$PC1, yvar=phylo_fem$PC2, factorvar = phylo_fem$habitat, labelvar = row.names(phylo_fem), edge.width = 0.5, title = 'b) Femelles', xlab = 'PC1 (61,8%)', ylab = 'PC2 (25,9%)') + scale_color_manual(values=mypalette, name='Habitat') + theme(legend.position = 'none', text = element_text(size=14))

######################### Signal phylogénétique Femelles #############################

Test.Kmult(mean_femelles, tree_morphos, iter=999)
# Signal phylogénétique Femelles moitié moindre que Males 

######################### MANNOVA phylogénétique Femelles #############################
b = ACP_mean_femelles$x[,1:7]
voltypr=as.factor(phylo_fem$habitat)
names(voltypr)<-dimnames(b)[[1]]
aov.phylo(b~voltypr, tree_morphos, nsim = 10, test="Wilks") 
summary(manova(b~voltypr))

############################## Figure phylomorphospaces ################################

grid_arrange_shared_legend(phylomorphomales,phylomorphofemelle, ncol=1, nrow=2, position='bottom')

################################# Dimorphisme sexuel ##################################

# On soustrait les valeurs moyennes de chaque espèce 
dimorphisme = mean_males-mean_femelles

# ACP sur dimorphisme
ACP_dimorphisme = prcomp(dimorphisme, scale. = T)
# Résumé de l'ACP
summary(ACP_dimorphisme)

### Graphes ACP 
# Explication des axes 
fviz_screeplot(ACP_dimorphisme, addlabels = T) 
# Représentation individus
fviz_pca_ind(ACP_dimorphisme, asp=1, pointshape=20, habillage=voltypr) + scale_color_discrete() + ggtitle('ACP dimorphisme sexuel')
fviz_pca_ind(ACP_dimorphisme, axes = c(2,3), asp=1, pointshape=20,  habillage=voltypr) + scale_color_discrete() + ggtitle('ACP dimorphisme sexuel') 

# Quantification du dimorphisme sexuel de chaque espèce en mesurant la norme entre les 2 sexes
norm_dimorphisme = sqrt((rowSums(dimorphisme^2)))
#barplot
toto<-barplot(norm_dimorphisme, col=mypalette[voltypr] ,axes = FALSE, axisnames = FALSE)
text(toto, par("usr")[3], labels = dimnames(dimorphisme)[[1]], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
axis(2)
title('Norme des vecteurs de dimorphisme')


########### Dimorphisme sexuel bis (femelles bleues cypris supprimées) ##############

# Suppression femelles cypris bleues
Femellesbis = spectres_femelles[,-44]
Femellesbis = Femellesbis[,-13]
Femellesbis = Femellesbis[,-114]
# Suppression des infos 
lignes = which(infos_femelles$ID=='MUSO'|infos_femelles$ID=='PBM68'|infos_femelles$ID=='PB13-0494')
infos_femelles_bis = infos_femelles[-lignes,]

# Calcul des spectres moyens 
spectres_mean_femelles_bis = aggspec(Femellesbis, by=infos_femelles_bis$espece, FUN = mean)
# Représentation graphique 
ggplot_rspec(spectres_mean_femelles_bis)

## ATTENTION ! la plupart des fonctions nécessitent variables en colonnes et individus en lignes
# C'est l'inverse des spectres (cf. modifications pour ACP)

# Transposition du dataframe
mean_femelles_bis = spec_to_table(spectres_mean_femelles_bis)

# On soustrait les valeurs moyennes de chaque espèce 
dimorphisme_bis = mean_males-mean_femelles_bis

# ACP sur dimorphisme
ACP_dimorphisme_bis = prcomp(dimorphisme_bis, scale. = T)
# Résumé de l'ACP
summary(ACP_dimorphisme_bis)

# Quantification du dimorphisme sexuel de chaque espèce en mesurant la norme entre les 2 sexes
norm_dimorphisme_bis = sqrt((rowSums(dimorphisme_bis^2)))
#barplot
toto_bis<-barplot(norm_dimorphisme_bis, col=mypalette[voltypr] ,axes = FALSE, axisnames = FALSE)
text(toto_bis, par("usr")[3], labels = dimnames(dimorphisme_bis)[[1]], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
axis(2)
title('Norme des vecteurs de dimorphisme')

# Signal phylogénétique
physignal(norm_dimorphisme_bis, tree_morphos, iter=999)

# Phylomorphospace
# Mise en couleur selon type de vol 
ACP_fem_jaune = prcomp(mean_femelles_bis, scale. = T)
summary(ACP_fem_jaune)
name.check(tree_morphos, ACP_fem_jaune$x[,1:2])
phylo_fem = as.data.frame(cbind(ACP_fem_jaune$x[,1:2]))

for (i in 1:nrow(phylo_fem)){
  canop = list(c('theseus', 'niepelti', 'amphitryon', 'telemachus', 'hercules', 'cisseis', 'hecuba', 'anaxibia', 'cypris', 'rhetenor'))
  if (row.names(phylo_fem[i,]) %in% canop[[1]] ){
    phylo_fem$habitat[i]='Canopée'
  }
  else {
    phylo_fem$habitat[i]='Sous-bois'
  }
}

phylomorphofemelle = ggphylomorpho(tree_morphos, phylo_fem, xvar = phylo_fem$PC1, yvar=phylo_fem$PC2, factorvar = phylo_fem$habitat, labelvar = row.names(phylo_fem), edge.width = 0.5, title = 'b) Femelles', xlab = 'PC1 (62,4%)', ylab = 'PC2 (26,5%)') + scale_color_manual(values=mypalette, name='Habitat') + theme(legend.position = 'none', text = element_text(size=14))

# Jointure mean Males & mean Femelles 
row.names(mean_males)=paste(row.names(mean_males), 'm')
row.names(mean_femelles_bis)=paste(row.names(mean_femelles_bis), 'f')
meanALL = rbind(mean_males, mean_femelles_bis) # .O = males, .1 = femelles
ACP_mean_all = prcomp(meanALL, scale. = T)
summary(ACP_mean_all)

infos_mean_all = c(rep('Male', 30), rep('Femelle', 30))

### Graphes ACP 
# Explication des axes 
fviz_screeplot(ACP_mean_all, addlabels = T) 
# Représentation variables 
#fviz_pca_var(ACP_mean_all, repel = T)
#fviz_pca_var(ACP_mean_all, axes = c(2,3), repel = T) 
# Représentation individus
fviz_pca_ind(ACP_mean_all, pointshape=20, pointsize=3, habillage = infos_mean_all, addEllipses = T, repel = T, labelsize=5, mean.point = FALSE, asp=1) + scale_color_discrete(name='Sexe') + ggtitle('') + theme(legend.position = c(0.95,0.1), legend.background = element_rect('white'), legend.text = element_text(size=12), text = element_text(size=14)) + labs(fill='Sexe') +xlim(-55, 55) + ylim (-55, 55)
fviz_pca_ind(ACP_mean_all, asp=1, axes = c(2,3), pointshape=20, habillage = infos_mean_all) + scale_color_discrete()

############# Figures rapport partie Macroévolution ################

#couleurs femelles avec cypris jaune
couleur_femelles = match.phylo.data(tree_morphos, mean_femelles_bis)$data
couleur_femelles = table_to_spec(couleur_femelles)

# Double phylogénie avec représentation approximative perception des spectres 
library(readxl)
both <- read_excel("~/Desktop/images.xlsx", sheet = "both") # Récupération path images
vol = c(rep('Canopy',10), rep('Understory', 20)) # type habitat selon tiplabels de tree_morpho
group = split(tree_morphos$tip.label, vol) # On sépare les tiplabels en fonction de l'habitat
tree2 = groupOTU(tree_morphos, group) # On cherche les branches de l'arbre reliant les 2 groupes C vs. U
# Arbre phylogénétique avec coloration branches + couleurs males 
t1=ggtree(tree2, ladderize=FALSE) + geom_tree(aes(color=group)) + geom_tippoint(aes(x=x+1), color=spec2rgb(as.rspec(couleur_males)), shape=19, size=5) + theme(legend.position = 'none') + scale_color_manual(values=c("#357AB7", "#82C46C"))
all = t1 %<+% both + xlim(0,20) # left-join path images + xlim pour libérer de la place à droite de la figure
# Ajout des couelurs des femelles, des tiplabels et des images 
all + geom_tippoint(aes(x=x+3), color=spec2rgb(as.rspec(couleur_femelles)), shape=19, size=5) + geom_tiplab(aes(image=imageURL_m), geom='image', size=0.035 , offset = 1.5) + geom_tiplab(aes(image=imageURL_f), geom='image', size=0.035 , offset = 3.5) + geom_tiplab(offset=5)


#######################################################################
############### TESTS DIMORPHISME BIJL ET AL 2020 ####################
#######################################################################

#### Librairies ####
library(ape)
library(tidytree)
library(devtools)
library(RRphylo)

# Fonction magnitude = dimorphisme 
## Van der Bijl et al. 2020 ##
magnitude <- function(x) {
  if (is.matrix(x)) sqrt(rowSums(x^2)) else sqrt(sum(x^2))
}

# Jeu de données 
jdd_dimorph = merge(mean_males, mean_femelles_bis, by='row.names', suffixes=c('males', 'femelles')) # x = males, y=femelles
rownames(jdd_dimorph)=jdd_dimorph$Row.names
jdd_dimorph=jdd_dimorph[,-1]
jdd_dimorph = as.matrix(jdd_dimorph)
# Vérif noms tree + jdd_dimorph
name.check(tree,jdd_dimorph)

# Estimation caractères aux noeuds en multivarié 
# Ridge regression method (Castiglione et al., 2018)
RR = RRphylo(tree, jdd_dimorph)

### Fonction permutations pour obtenir la distribution nulle (H0) 
## Van der Bijl et al. 2020 ##
fit_permuted_RR <- function(m, tree) {
  require(RRphylo)
  # randomly choose whether to keep or swap the male/female label
  r <- sample(c(TRUE, FALSE), nrow(m), replace = TRUE)
  m_swapped <- m[, c(402:802, 1:401)]
  colnames(m_swapped) <- colnames(m)
  
  # create a permuted matrix and reorder
  m_perm <- rbind(m[!r, ], m_swapped[r, ])
  m_perm <- m_perm[tree$tip.label, ]
  RRphylo(tree, m_perm, clus = 0)
}

library(future.apply)
# permuations 
#plan(multisession, workers = 16) 
#perms <- future_replicate(1e3, fit_permuted_RR(jdd_dimorph, tree), FALSE, FALSE)
# Sauver en .rd car long à lancer à chaque fois 
#write_rds(perms, '~/Desktop/evolrate_RR_perms.rds')
perm <- read_rds('~/Desktop/evolrate_RR_perms.rds')
# ratio of evolutionary rates in the permutations:
p_diffs <- map_dbl(perm, ~mean(magnitude(.x$multiple.rates[,1:401]) / mean(magnitude(.x$multiple.rates[,402:802]))))
# observed ratio of evolutionary rates:
o_diff <- mean(magnitude(RR$multiple.rates[,1:401])) / mean(magnitude(RR$multiple.rates[,402:802]))
# p-value:
pmin(mean(p_diffs >= o_diff), mean(p_diffs <= o_diff)) * 2
# perm<=obs n'existe pas => tous perm > obs donc p=O.OO1

## Summarise as figure. S2 in the paper.
## Van der Bijl et al. 2020 ##
ggplot() + 
  geom_vline(aes(xintercept = 1), size = 1, color = 'grey50') +
  geom_density(aes(x, fill = 'simulé'), data.frame(x = p_diffs), alpha = 0.6, col = NA) + 
  geom_vline(aes(xintercept = x, color = 'observée'), data.frame(x = o_diff), size = 1) +
  scale_fill_manual(values = 'grey20') +
  scale_x_log10(
    limits = c(7.5/10, 10/7.5), breaks = c(4/5, 10/9, 1, 9/10, 5/4), 
    labels = expression(frac(4, 5), frac(10, 9), 1, frac(9, 10), frac(5, 4))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = 'firebrick') +
  labs(
    x = expression(frac('taux d\'évolution male', 'taux d\'évolution femelle')),
    y = 'densité',
    lty = NULL, color = NULL, fill = NULL,
    caption = paste('p =', pmin(mean(o_diff <= p_diffs), mean(o_diff >= p_diffs)) * 2)
  ) +
  theme_minimal() + 
  theme(legend.position = 'top', legend.margin = margin(0, 0, 0, 0))
ggsave('~/Desktop/fig1.png', dpi = 400, width = 5, height = 2.5)

# Get the estimated ancestral states, combined with tip values
aces_m <- rbind(RR$aces[, 1:401], jdd_dimorph[, 1:401])
aces_f <- rbind(RR$aces[, 402:802], jdd_dimorph[, 402:802])

# dichromatism:
dE <- sqrt(rowSums((aces_f - aces_m)^2))
# evolutionary rates:
rF <- magnitude(RR$multiple.rates[, 402:802])
rM <- magnitude(RR$multiple.rates[, 1:401])

# Function to obtain the correlation between dichromatism and evolutionary rates.
extract_cors <- function(jdd_dimorph, RR, plot = FALSE) {
  aces_m <- rbind(RR$aces[, 1:401], jdd_dimorph[, 1:401])
  aces_f <- rbind(RR$aces[, 402:802], jdd_dimorph[, 402:802])
  dE <- magnitude(aces_f - aces_m)
  rF <- magnitude(RR$multiple.rates[, 402:802])
  rM <- magnitude(RR$multiple.rates[, 1:401])
  
  return(c(
    ratio  = coef(lm(log(rM / rF) ~ dE))[2],
    ratio_intercept = coef(lm(log(rM / rF) ~ dE))[1]
  ))
}

# observed slope
O_ratio <- extract_cors(jdd_dimorph, RR)[2]
# expected slopes
E_ratio <- map_dbl(perm, ~extract_cors(jdd_dimorph, .)[2])
# p-value
p_ratio <- pmin(mean(O_ratio >= E_ratio), mean(O_ratio <= E_ratio)) * 2
p_ratio
if (p_ratio == 0) p_ratio <- 0.001


# Make figure 3:
p1d <- data_frame(dE, rF = rF, rM = rM)

p1 <- ggplot(p1d, aes(dE, rM / rF)) +
  geom_point(size = 1, shape = 21, alpha = 0.5, col = 1, fill = 1) + 
  geom_hline(yintercept = 1, color = 'grey60') +
  geom_smooth(method = 'lm', se = FALSE, show.legend = FALSE, col = 'firebrick') +
  scale_y_log10() +
  labs(x = 'dichromatism (||D||)', y = expression(frac('male rate', 'female rate'))) +
  theme_minimal()

xmax <- 2500
p1
p1_plus <- p1 + 
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    map_dfr(perm, function(.) {
      x <- extract_cors(jdd_dimorph, .)
      data.frame(x = 0, y = exp(x[2]), xend = xmax, yend = exp(xmax * x[1] + x[2]))
    } ), 
    size = 0.05, alpha = 0.05
  ) +
  geom_smooth(method = 'lm', se = FALSE, show.legend = FALSE, col = 'firebrick')
ggsave('~/Desktop/figure3.png', p1_plus, height = 4, width = 4, dpi = 400)
ggsave('~/Desktop/figure3.pdf', p1_plus, height = 4, width = 4)

scalar_projection <- function(a, b) { 
  sca_proj <- function(A, B) (A %*% B) / sqrt(sum(B^2))
  if (is.vector(a) & is.vector(b)) return(sca_proj(a, b))
  map_dbl(1:nrow(a), ~sca_proj(a[., ], b[., ]))
}

### Prep -------------------------------------------------------------------------------------------

# The predicted Euclidean distance between male and female phenotypes (i.e. dichromatism)
E <- sqrt(rowSums((aces_m - aces_f)^2))
# Evolutionary vectors in 3d color space
vF <- RR$multiple.rates[, 402:802]
vM <- RR$multiple.rates[, 1:401]
# The vectors from female to male phenotype, i.e. direction and size of dichromatism
vE <- aces_m - aces_f
# The evolutionary rate of dichromatism itself
rE <- sqrt(rowSums(vE ^ 2))
# Evolutionary rates for males and females, i.e. the vector magnitude
rF <- magnitude(RR$multiple.rates[, 402:802])
rM <- magnitude(RR$multiple.rates[, 1:401])
# Effective rates for males and females, i.e. the size of the projection on the the dichromatic direction
rF_eff <- scalar_projection(-vF, vE)
rM_eff <- scalar_projection(vM, vE)
# Vector of change in dichromatism
vdE <- vM - vF
rvdE <- sqrt(rowSums((vdE)^2))
rE_eff <- scalar_projection(vdE, vE)

test_plot = as.data.frame(cbind(rM,rF))
test_plot$habillage = rM/rF
ggplot(test_plot, aes(x=rF, y=rM, color=ifelse(habillage>1, T, 'green'))) + 
  geom_point() + theme(legend.position = "none") + 
  geom_segment(
  aes(x = 0, y = 0, xend = 150, yend = 150, color='black')) + 
  coord_equal()

### Fit and analyze permuted -----------------------------------------------------------------------

extract_cors <- function(rr) {
  aces_f <- rbind(rr$aces[, 402:802], jdd_dimorph[, 402:802])
  aces_m <- rbind(rr$aces[, 1:401], jdd_dimorph[, 1:401])
  dE <- magnitude(aces_f - aces_m)
  # Estimated evolutionary changes per axis (3d vectors)
  vF <- rr$multiple.rates[, 402:802]
  vM <- rr$multiple.rates[, 1:401]
  # Estimated ancestral state of dichromatism
  vE <- aces_m - aces_f
  # The evolutionary rate of dichromatism
  rE <- sqrt(rowSums(vE ^ 2))
  # Evolutionary rates for males and females, i.e. the vector magnitude
  rF <- magnitude(rr$multiple.rates[, 402:802])
  rM <- magnitude(rr$multiple.rates[, 1:401])
  # Effective rates for males and females, i.e. the size of the projection onto the dichromatic direction
  rF_eff <- scalar_projection(-vF, vE)
  rM_eff <- scalar_projection(vM, vE)
  # Vector of change in dichromatism
  vdE <- vM - vF
  rF_eff2 <- scalar_projection(-vF, vdE)
  rM_eff2 <- scalar_projection(vM, vdE)
  # Effictive rate of change in dichromatism, i.e. change in dichromatism projected onto existing dichromatism
  rE_eff <- scalar_projection(vdE, vE)
  
  m1 <- lm(rM_eff ~ rE_eff)
  m2 <- lm(rF_eff ~ rE_eff)
  
  return(setNames(
    c(coef(m1), coef(m2)),
    c('male_intercept', 'male_slope', 'female_intercept', 'female_slope')
  ))
}

extract_cors(RR)

library(furrr)

# observed difference in slopes
O <- diff(extract_cors(RR)[c(2, 4)])
# expected difference in slopes
plan(multisession, workers = 16)
E <- future_map_dbl(perm, ~diff(extract_cors(.)[c(2, 4)]), .progress = TRUE)
plan(sequential)
# p-value
p <- pmin(mean(O >= E), mean(O <= E)) * 2
if (p == 0) p <- 0.001

# Make figure 4
p1 <- ggplot(mapping = aes(rE_eff, rM_eff)) +
  geom_point(alpha = 0.3, shape = 21, fill = 1) +
  theme_minimal()
xmin <- min(rE_eff); xmax <- max(rE_eff)
p1_plus <- p1 + 
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    map_dfr(perm, function(.) {
      x <- extract_cors(.)
      data.frame(x = xmin, y = xmin * x[2] + x[1], xend = xmax, yend = xmax * x[2] + x[1])
    } ), 
    size = 0.05, alpha = 0.025
  ) +
  geom_smooth(method = 'lm', se = FALSE, col = 'firebrick') +
  labs(x = 'Évolution du dimorphisme', y = 'Part d\'évolution du dimorsphisme\ndue aux mâles')

p2 <- p1 + aes(rE_eff, rF_eff)
p2_plus <- p2 + 
  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    map_dfr(perm, function(.) {
      x <- extract_cors(.)
      data.frame(x = xmin, y = xmin * x[4] + x[3], xend = xmax, yend = xmax * x[4] + x[3])
    } ), 
    size = 0.05, alpha = 0.025
  ) +
  geom_smooth(method = 'lm', se = FALSE, col = 'firebrick') +
  labs(x = 'Évolution du dimorphisme', y = 'Part d\'évolution du dimorsphisme\ndue aux femelles')

cowplot::plot_grid(p1_plus, p2_plus, labels = 'auto')
ggsave('~/Desktop/figure4.pdf', height = 4, width = 7)

###################### Jeu de données vol #######################

setwd('~/Desktop/Vol')

# Importation du jeu de données vol 
flight_parameters = read.csv('~/Desktop/vol/0_datavol.csv', header=T, sep=";", na.strings = NA)
supprimer = which(flight_parameters$specimen_ID==79|flight_parameters$specimen_ID==92)
flight_parameters = flight_parameters[-supprimer,]

# Importation des spectres 
spectra_vol = getspec("~/Desktop/Vol", ext = "RFL8", lim = c(300,700), subdir = F)
# Localisation des fichiers 
specfiles_vol = list.files(path='.', pattern='*\\.RFL8',full.names = TRUE, recursive = TRUE)
# Récupération des infos 
infos_vol = do.call(rbind, lapply(specfiles_vol, recup_infos))
# smooth
spectra_vol_sm = procspec(spectra_vol, opt = 'smooth', fixneg = 'zero', span = 0.2)

spectres_vol = spec_to_table(spectra_vol_sm)
row.names(spectres_vol)=infos_vol$ID
spectres_vol=as.matrix(spectres_vol)

# Traitement NA dans paramètres de vol 
library(missMDA)

flight_parameters = cbind(infos_vol$ID, infos_vol$espece, flight_parameters$mean_flight_height, flight_parameters$mean_curvature, flight_parameters$sinuosity, flight_parameters$sink_speed, flight_parameters$mean_acceleration, flight_parameters$glide_ratio, flight_parameters$Micro_wb_frequency, flight_parameters$relative_advance_ratio, flight_parameters$mean_ascent_angle, flight_parameters$mean_glide_duration, flight_parameters$glide_angle)
colnames(flight_parameters)=c('ID', 'sp', 'mean_flight_height', 'mean_curvature', 'sinuosity', 'sink_speed', 'mean_acceleration', 'glide_ratio', 'Micro_wb_frequency', 'relative_advance_ratio', 'mean_ascent_angle', 'mean_glide_duration', 'glide_angle')
flight_parameters=as.data.frame(flight_parameters)
row.names(flight_parameters)=flight_parameters$ID
flight_parameters=flight_parameters[,-1]
flight_parameters=flight_parameters[,-1]

sapply(flight_parameters, class)
i = c(1:11)
flight_parameters[ , i] <- apply(flight_parameters[ , i], 2, function(x) as.numeric(as.character(x)))

#### remplace les NA avec des estimations
imputed_data <- imputePCA(flight_parameters) # re-estimate the NAs values
#### les deux lignes suivantes c'est du formatage
imputed_data <- imputed_data$completeObs
imputed_data <- as.data.frame(imputed_data)

# Suppression n=1 + femelle 
spectra_vol_sm = spectra_vol_sm[,-c(37,52)]
imputed_data = imputed_data[-c(36,51),]
infos_vol=infos_vol[-c(36,51),]
spectres_vol = spectres_vol[-c(36,51),]
row.names(spectres_vol)=row.names(imputed_data)

indv_param_col_vol = summary(spectra_vol_sm, subset=c('H1','B1', 'S2', 'S3'), wlmin=325)
row.names(indv_param_col_vol)=infos_vol$ID

### Corrplot ####
library(ggcorrplot)
matrice_indv = cbind(indv_param_col_vol, imputed_data)
M_indv = cor(matrice_indv)
pmat_indv = cor_pmat(M_indv)
for(i in 1:nrow(M_indv)){
  for(j in 1:nrow(M_indv)){
    if(pmat_indv[i,j] > 0.05){
      M_indv[i,j] <- NA
    }
  }
}

ggcorrplot(M_indv, method = 'circle')

library(RColorBrewer)
mypalette=brewer.pal(10, "Spectral")


b = prcomp(spectres_vol, scale. = T)
summary(b)

plot(prcomp(imputed_data)$x, pch=20, cex=2, col= mypalette[factor(infos_vol$espece)])
text(prcomp(imputed_data)$x, labels=infos_vol$espece, pos=3)

PLS_non_phylo = two.b.pls(b$x[,1:9], imputed_data, iter=999)#PLS non phylogenetiques
PLS_non_phylo
plot(PLS_non_phylo, label=infos_vol$espece) 


# moyennes d'espèces
mean_spectre_vol = aggspec(spectra_vol_sm, by=infos_vol$espece, FUN=mean)
# Paramètres de couleur 
param_col_vol = summary(mean_spectre_vol, subset=c('H1','B1', 'S2', 'S3'), wlmin=325)
# Transposition du dataframe
mean_spectre_vol = spec_to_table(mean_spectre_vol)

c=prcomp(mean_spectre_vol, scale. = T)
summary(c)
mean_param_vol = meangr(imputed_data, infos_vol$espece)

#si jamais mismatch tree data: on coupe l'arbre
NAs<-name.check(tree_morphos,mean_param_vol)$tree_not_data
tree2<-drop.tip(tree_morphos,(NAs))
name.check(tree2, mean_param_vol)

PLS_non_phylo = two.b.pls(c$x[,1:7], mean_param_vol, iter=999)#PLS non phylogenetiques
PLS_non_phylo
plot(PLS_non_phylo, label=row.names(mean_param_vol)) 

PLS_phylo = phylo.integration(c$x[,1:7], mean_param_vol, tree2, iter=999)
PLS_phylo
plot(PLS_phylo, label = row.names(mean_param_vol))


#### Corrplot ####

matrice = cbind(param_col_vol, mean_param_vol)
M = cor(matrice)
pmat = cor_pmat(M)
for(i in 1:nrow(M)){
   for(j in 1:nrow(M)){
     if(pmat[i,j] > 0.05){
       M[i,j] <- NA
     }
   }
 }
  
ggcorrplot(M, method = 'circle')


plot(param_col_vol$B1~mean_param_vol[,1])
text(param_col_vol$B1~mean_param_vol[,1], labels=row.names(param_col_vol), pos=2)

plot(param_col_vol$B1~mean_param_vol[,11])
text(param_col_vol$B1~mean_param_vol[,11], labels=row.names(param_col_vol), pos=2)

plot(param_col_vol$S3~mean_param_vol[,2])
text(param_col_vol$S3~mean_param_vol[,2], labels=row.names(param_col_vol), pos=2)

#PLS non phylo
PLS_non_phylo = two.b.pls(param_col_vol, mean_param_vol, iter=999)#PLS non phylogenetiques
PLS_non_phylo
plot(PLS_non_phylo, label=row.names(mean_param_vol)) 

#PLS phylo
PLS_phylo = phylo.integration(param_col_vol, mean_param_vol, tree2, iter=999)
PLS_phylo
plot(PLS_phylo, label = row.names(mean_param_vol))

############## Sous jeu de données vol ############## 

### ACHILLES
achilles_vol = which(infos_vol$espece=='achilles')
achilles_param_vol = imputed_data[achilles_vol,]
achilles_param_col = indv_param_col_vol[achilles_vol,]

matrice_achilles = cbind(achilles_param_col, achilles_param_vol)
M_achilles = cor(matrice_achilles)
pmat_achilles = cor_pmat(M_achilles)
for(i in 1:nrow(M_achilles)){
  for(j in 1:nrow(M_achilles)){
    if(pmat_achilles[i,j] > 0.05){
      M_achilles[i,j] <- NA
    }
  }
}

ggcorrplot(M_achilles, method = 'circle')

PLS_achilles = two.b.pls(achilles_param_col, achilles_param_vol, iter=999)#PLS non phylogenetiques
PLS_achilles
plot(PLS_achilles)

### HELENOR
helenor_vol = which(infos_vol$espece=='helenor')
helenor_param_vol = imputed_data[helenor_vol,]
helenor_param_col = indv_param_col_vol[helenor_vol,]

matrice_helenor = cbind(helenor_param_col, helenor_param_vol)
M_helenor = cor(matrice_helenor)
pmat_helenor = cor_pmat(M_helenor)
for(i in 1:nrow(M_helenor)){
  for(j in 1:nrow(M_helenor)){
    if(pmat_helenor[i,j] > 0.05){
      M_helenor[i,j] <- NA
    }
  }
}

ggcorrplot(M_helenor, method = 'circle')

PLS_helenor = two.b.pls(helenor_param_col, helenor_param_vol, iter=999)#PLS non phylogenetiques
PLS_helenor
plot(PLS_helenor)

### GODARTII
godartii_vol = which(infos_vol$espece=='godartii')
godartii_param_vol = imputed_data[godartii_vol,]
godartii_param_col = indv_param_col_vol[godartii_vol,]

matrice_godartii = cbind(godartii_param_col, godartii_param_vol)
M_godartii = cor(matrice_godartii)
pmat_godartii = cor_pmat(M_godartii)
for(i in 1:nrow(M_godartii)){
  for(j in 1:nrow(M_godartii)){
    if(pmat_godartii[i,j] > 0.05){
      M_godartii[i,j] <- NA
    }
  }
}

ggcorrplot(M_godartii, method = 'circle')

PLS_godartii = two.b.pls(godartii_param_col, godartii_param_vol, iter=999)#PLS non phylogenetiques
PLS_godartii
plot(PLS_godartii)


### THESEUS
theseus_vol = which(infos_vol$espece=='theseus')
theseus_param_vol = imputed_data[theseus_vol,]
theseus_param_col = indv_param_col_vol[theseus_vol,]

matrice_theseus = cbind(theseus_param_col, theseus_param_vol)
M_theseus = cor(matrice_theseus)
pmat_theseus = cor_pmat(M_theseus)
for(i in 1:nrow(M_theseus)){
  for(j in 1:nrow(M_theseus)){
    if(pmat_theseus[i,j] > 0.05){
      M_theseus[i,j] <- NA
    }
  }
}

ggcorrplot(M_theseus, method = 'circle')

PLS_theseus = two.b.pls(theseus_param_col, theseus_param_vol, iter=999)#PLS non phylogenetiques
PLS_theseus
plot(PLS_theseus)

### RHETENOR
rhetenor_vol = which(infos_vol$espece=='rhetenor')
rhetenor_param_vol = imputed_data[rhetenor_vol,]
rhetenor_param_col = indv_param_col_vol[rhetenor_vol,]

matrice_rhetenor = cbind(rhetenor_param_col, rhetenor_param_vol)
M_rhetenor = cor(matrice_rhetenor)
pmat_rhetenor = cor_pmat(M_rhetenor)
for(i in 1:nrow(M_rhetenor)){
  for(j in 1:nrow(M_rhetenor)){
    if(pmat_rhetenor[i,j] > 0.05){
      M_rhetenor[i,j] <- NA
    }
  }
}

ggcorrplot(M_rhetenor, method = 'circle')

PLS_rhetenor = two.b.pls(rhetenor_param_col, rhetenor_param_vol, iter=999)#PLS non phylogenetiques
PLS_rhetenor
plot(PLS_rhetenor)

### MENELAUS
menelaus_vol = which(infos_vol$espece=='menelaus')
menelaus_param_vol = imputed_data[menelaus_vol,]
menelaus_param_col = indv_param_col_vol[menelaus_vol,]

matrice_menelaus = cbind(menelaus_param_col, menelaus_param_vol)
M_menelaus = cor(matrice_menelaus)
pmat_menelaus = cor_pmat(M_menelaus)
for(i in 1:nrow(M_menelaus)){
  for(j in 1:nrow(M_menelaus)){
    if(pmat_menelaus[i,j] > 0.05){
      M_menelaus[i,j] <- NA
    }
  }
}

ggcorrplot(M_menelaus, method = 'circle')

PLS_menelaus = two.b.pls(menelaus_param_col, menelaus_param_vol, iter=999)#PLS non phylogenetiques
PLS_menelaus
plot(PLS_menelaus)

### DEIDAMIA
deidamia_vol = which(infos_vol$espece=='deidamia')
deidamia_param_vol = imputed_data[deidamia_vol,]
deidamia_param_col = indv_param_col_vol[deidamia_vol,]

matrice_deidamia = cbind(deidamia_param_col, deidamia_param_vol)
M_deidamia = cor(matrice_deidamia)
pmat_deidamia = cor_pmat(M_deidamia)
for(i in 1:nrow(M_deidamia)){
  for(j in 1:nrow(M_deidamia)){
    if(pmat_deidamia[i,j] > 0.05){
      M_deidamia[i,j] <- NA
    }
  }
}

ggcorrplot(M_deidamia, method = 'circle')

PLS_deidamia = two.b.pls(deidamia_param_col, deidamia_param_vol, iter=999)#PLS non phylogenetiques
PLS_deidamia
plot(PLS_deidamia)

### CISSEIS
cisseis_vol = which(infos_vol$espece=='cisseis')
cisseis_param_vol = imputed_data[cisseis_vol,]
cisseis_param_col = indv_param_col_vol[cisseis_vol,]

matrice_cisseis = cbind(cisseis_param_col, cisseis_param_vol)
M_cisseis = cor(matrice_cisseis)
pmat_cisseis = cor_pmat(M_cisseis)
for(i in 1:nrow(M_cisseis)){
  for(j in 1:nrow(M_cisseis)){
    if(pmat_cisseis[i,j] > 0.05){
      M_cisseis[i,j] <- NA
    }
  }
}

ggcorrplot(M_cisseis, method = 'circle')

PLS_cisseis = two.b.pls(cisseis_param_col, cisseis_param_vol, iter=999)#PLS non phylogenetiques
PLS_cisseis
plot(PLS_cisseis)

### SULKOWSKYI => n=2
### MARCUS => n=2

###################### Jeu de données sympatrie #######################

setwd('~/Desktop/Sympatry')

# Localisation des fichiers 
specfiles = list.files(path='.', pattern='*\\.RFL8',full.names = TRUE, recursive = TRUE)
# Récupération des infos 
infos = do.call(rbind, lapply(specfiles, recup_infos))
# Importation des spectres de 300 à 700 nm
spectra_sympat = getspec("~/Desktop/Sympatry", ext = "RFL8", lim = c(300,700), subdir = F)
# Que les bleus 
bleu_sympat = subset(spectra_sympat, 'bleu')

# Smooth de tous les spectres avec le même coefficient 
spectra_sym_sm = procspec(bleu_sympat, opt = 'smooth', fixneg = 'zero', span = 0.2)

#Males des 3 espèces 
symp_males = subset(spectra_sym_sm, 'Male')
# Localisation des fichiers 
specfiles = list.files(path='.', pattern='*\\Male_bleu.RFL8',full.names = TRUE, recursive = TRUE)
# Récupération des infos 
symp_bleu = do.call(rbind, lapply(specfiles, recup_infos))
symp_bleu$Loc = c(rep('Pérou',31), rep('Guyanne',20))
both = paste(symp_bleu$espece, symp_bleu$Loc)

AchillesM = subset(symp_males, 'achilles')
HelenorM = subset(symp_males, 'helenor')
DeidamiaM = subset(symp_males, 'deidamia')

# Achilles
infos_achilles_M = c(rep('Pérou',11), rep('Guyanne', 9))
aggplot(AchillesM, by=infos_achilles_M, legend = T, main = 'Achilles')

# Helenor
infos_helenor_M = c(rep('Pérou',10), rep('Guyanne', 6))
aggplot(HelenorM, by=infos_helenor_M, legend = T, main = 'Helenor')

# Deidamia
infos_deidamia_M = c(rep('Pérou',10), rep('Guyanne', 5))
aggplot(DeidamiaM, by=infos_deidamia_M, legend = T, main = 'Deidamia')

### ACP 
ACP = prcomp(spec_to_table(symp_males), scale. = T)
# Résumé de l'ACP
summary(ACP)

# Graphes ACP 

fviz_screeplot(ACP, addlabels = T)
fviz_pca_ind(ACP, asp=1, labels = F, pointsize=0) + geom_point(aes(color=symp_bleu$espece, shape=symp_bleu$Loc), size=2) + scale_shape_manual(values=c(16,21))
fviz_pca_ind(ACP, axes = c(1,3), asp=1, labels = F, pointsize=0) + geom_point(aes(color=symp_bleu$Loc, shape=symp_bleu$espece), size=2) + scale_shape_manual(values=c(15, 16, 17))

param_males = summary(symp_males, subset=c('H1','B1', 'S2', 'S3'), wlmin=325)
hele = summary(HelenorM, subset=c('H1','B1', 'S2', 'S3'), wlmin=325)
deid = summary(DeidamiaM, subset=c('H1','B1', 'S2', 'S3'), wlmin=325)
achill = summary(AchillesM, subset=c('H1','B1', 'S2', 'S3'), wlmin=325)

# MANOVA 
summary(manova(ACP$x[,1:7]~symp_bleu$espece))
summary(manova(ACP$x[,1:7]~symp_bleu$Loc))
summary(manova(ACP$x[,1:7]~symp_bleu$Loc*symp_bleu$espece))

# Boxplots paramètres 
param_males=cbind(param_males, symp_bleu$espece, symp_bleu$Loc)
p1 = ggplot(param_males, aes(x=symp_bleu$espece, y=H1, fill=symp_bleu$Loc)) + geom_boxplot() + theme(legend.position = 'none') + xlab('') + ylab('Teinte (H1)') + theme_classic() + theme(axis.text = element_text(size=12), axis.title = element_text(size=14))
p2 = ggplot(param_males, aes(x=symp_bleu$espece, y=B1, fill=symp_bleu$Loc)) + geom_boxplot() + theme(legend.position = 'none') + xlab('') + ylab('Luminance (B1)')+ theme_classic() + theme(axis.text = element_text(size=12), axis.title = element_text(size=14))
p3 = ggplot(param_males, aes(x=symp_bleu$espece, y=S2, fill=symp_bleu$Loc)) + geom_boxplot() + theme(legend.position = 'none') + xlab('') + ylab('Saturation (S2)')+ theme_classic() + theme(axis.text = element_text(size=12), axis.title = element_text(size=14))
p4 = ggplot(param_males, aes(x=symp_bleu$espece, y=S3, fill=symp_bleu$Loc)) + geom_boxplot() + theme(legend.position = 'none') + xlab('') + ylab('Chroma (S3)')+ theme_classic() + theme(axis.text = element_text(size=12), axis.title = element_text(size=14))
grid_arrange_shared_legend(p1,p2, p3, p4, ncol=2, nrow=2, position='bottom') 

# ANOVAs
summary(aov(param_males$H1~symp_bleu$espece)) #NS
# vérif conditions d'application 
a=aov(param_males$H1~symp_bleu$espece)
qqnorm(a$residuals) # Normalité des résidus
bartlett.test(a$residuals, symp_bleu$espece) # Homoscédasticité

summary(aov(param_males$H1~symp_bleu$Loc)) #***
summary(aov(param_males$H1~symp_bleu$Loc*symp_bleu$espece)) #NS
summary(aov(param_males$B1~symp_bleu$espece)) #***
summary(aov(param_males$B1~symp_bleu$Loc)) #NS
summary(aov(param_males$B1~symp_bleu$Loc*symp_bleu$espece)) #**
summary(aov(param_males$S2~symp_bleu$espece)) #***
summary(aov(param_males$S2~symp_bleu$Loc)) #NS
summary(aov(param_males$S2~symp_bleu$Loc*symp_bleu$espece)) #***
summary(aov(param_males$S3~symp_bleu$espece)) #***
summary(aov(param_males$S3~symp_bleu$Loc)) #*
summary(aov(param_males$S3~symp_bleu$Loc*symp_bleu$espece)) #*


#### Tests vision models ####

#### Fonctions d'affichage des graphes ####
graphe_chro = function(x, y){
  tableau=as.data.frame(x)
  #lignes = c(9,7,14,4,2,11,1,15,10)
  #groupes = c(rep('Pérou',3), rep('Guyane', 3), rep('Intraspécifique',3))
  #associations = c(rep(c('achilles-helenor', 'achilles-deidamia', 'helenor-deidamia'),2), c('achilles-achilles', 'helenor-helenor', 'deidamia-deidamia'))
  #tableau = x[lignes,]
  #tableau = as.data.frame(cbind(tableau, groupes, associations))
  #i = c(1:6)
  #tableau[ , i] <- apply(tableau[ , i], 2, function(x) as.numeric(as.character(x)))
  #row.names(tableau)=c(1:9)
  #ymaxdS = as.numeric(max(tableau[,1])) + as.numeric(max(tableau[,1]))/3
  ggplot(tableau) + 
    stat_sum(aes(x=row.names(tableau), y=dS.mean, ymin=dS.lwr, ymax=dS.upr), geom = "pointrange", shape=20, size=1) + 
    scale_x_discrete(labels=row.names(tableau)) +
    geom_abline(intercept=1,slope=0,linetype='dashed')+
    ylab('Contraste chromatique (en JND)')+
    xlab('')+
    ggtitle(y)+
    theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = 'right', plot.title = element_text(size=13, face='bold', hjust = 0.5))
}

graphe_achro = function(x,y){
  tableau=as.data.frame(x)
  #lignes = c(9,7,14,4,2,11,1,15,10)
  #groupes = c(rep('Pérou',3), rep('Guyane', 3), rep('Intraspécifique',3))
  #associations = c(rep(c('achilles-helenor', 'achilles-deidamia', 'helenor-deidamia'),2), c('achilles-achilles', 'helenor-helenor', 'deidamia-deidamia'))
  #tableau = x[lignes,]
  #tableau = as.data.frame(cbind(tableau, groupes, associations))
 # i = c(1:6)
  #tableau[ , i] <- apply(tableau[ , i], 2, function(x) as.numeric(as.character(x)))
  #row.names(tableau)=c(1:9)
 # ymaxdL = as.numeric(max(tableau[,4])) + as.numeric(max(tableau[,4]))/3
  ggplot(tableau) + 
    stat_sum(aes(x=row.names(tableau), y=dL.mean, ymin=dL.lwr, ymax=dL.upr), geom = "pointrange", shape=20, size=1) + 
    scale_x_discrete(labels=row.names(tableau)) +
    geom_abline(intercept=1,slope=0,linetype='dashed')+
    ylab('Contraste achromatique (en JND)')+
    xlab('')+
    ggtitle(y)+
    theme(axis.text.x = element_text(angle=60, hjust=1), legend.position = 'right', plot.title = element_text(size=13, face='bold', hjust = 0.5))
}

# Modèle oiseaux UV-sensibles
visualmodUV = vismodel(symp_males, visual = 'avg.uv', achromatic='bt.dc', illum = 'forestshade')
espacecolUV = colspace(visualmodUV, space='tcs')
entrespetpop = bootcoldist(visualmodUV, n=c(0.37,0.7,0.99,1), weber=0.05, weber.achro=0.05, by=both)
entrespetpop
graphe_achro(entrespetpop, 'Oiseau UV sensible')
graphe_chro(entrespetpop, 'Oiseau UV sensible')

# Modèle oiseaux V-sensibles
visualmodV = vismodel(symp_males, visual = 'avg.v', achromatic='ch.dc', illum = 'forestshade') #illum = 'forestshade'
espacecolV = colspace(visualmodV, space='tcs')
entrespetpopV = bootcoldist(visualmodV, n=c(0.25,0.5,1,1),weber=0.05,weber.achro=0.5, by=both)
entrespetpopV
graphe_achro(entrespetpopV, 'Oiseau Violet sensible')
graphe_chro(entrespetpopV, 'Oiseau Violet sensible')

#### Heliconius visual system ####
H_erato_sens = sensmodel(c(355, 390, 470, 555))
H_erato_vis = vismodel(symp_males, visual = H_erato_sens, relative = F, achromatic ='all', illum = 'forestshade')
H_erato_colspace = colspace(H_erato_vis, space='tcs')
H_erato_dist_e = bootcoldist(H_erato_vis, by=both, n=c(0.09, 0.07, 0.17, 1), weber=0.05, weber.achro=0.05) #Finkbeiner et al., 2017
H_erato_dist_e
graphe_achro(H_erato_dist_e, 'Heliconius erato femelle')
graphe_chro(H_erato_dist_e, 'Heliconius erato femelle')

#### E. atala (Lycaenidae) visual system ####
## D'après Liénard et al., 2021 ##

E_atala_sens = sensmodel(c(362, 441, 495, 564))
E_atala_vis = vismodel(symp_males, visual = E_atala_sens, relative = F, achromatic = 'all', illum = 'forestshade')
E_atala_colspace = colspace(E_atala_vis, space='tcs')
E_atala_dist_e = bootcoldist(E_atala_vis, by=both, n=c(0.06, 0.16, 0.03, 0.75), weber=0.05, weber.achro=0.05) #Finkbeiner et al., 2017
E_atala_dist_e
graphe_achro(E_atala_dist_e, 'Eumaeus atala')
graphe_chro(E_atala_dist_e, 'Eumaeus atala')

#### Chromatic contrast Homo sapiens ####
visualmodH = vismodel(symp_males, visual = 'cie2', illum = 'forestshade')
espacecolH = colspace(visualmodH, space='tri')
entrespetpopH = bootcoldist(espacecolH, by=both)
entrespetpopH

#### Fonction pour graphe tout ensemble ###
graphe = function(w, wname, x, xname, y, yname, z, zname, h, hname){
  lignes = c(9,7,14,4,2,11,1,15,10)
  groupes = c(rep('Pérou',3), rep('Guyane', 3), rep('Interlocalités',3))
  associations = c(rep(c('achilles-helenor', 'achilles-deidamia', 'helenor-deidamia'),2), c('achilles-achilles', 'helenor-helenor', 'deidamia-deidamia'))
  i = c(1:6)
  tableau_w = w[lignes,]
  tableau_w = as.data.frame(cbind(tableau_w, groupes, associations, Observateur = as.character(rep(wname,length(lignes)))))
  tableau_x = x[lignes,]
  tableau_x = as.data.frame(cbind(tableau_x, groupes, associations, Observateur = as.character(rep(xname,length(lignes)))))
  tableau_y = y[lignes,]
  tableau_y = as.data.frame(cbind(tableau_y, groupes, associations, Observateur = as.character(rep(yname,length(lignes)))))
  tableau_z = z[lignes,]
  tableau_z = as.data.frame(cbind(tableau_z, groupes, associations, Observateur = as.character(rep(zname,length(lignes)))))
  tableau_h = h[lignes,]
  tableau_h = as.data.frame(cbind(tableau_h, groupes, associations, Observateur = as.character(rep(hname,length(lignes)))))
  contrastes = as.data.frame(rbind(tableau_w, tableau_x, tableau_y, tableau_z, tableau_h))
  contrastes = as.data.frame(contrastes[c(which(contrastes$associations=='achilles-helenor' & contrastes$groupes=='Pérou'), 
                            which(contrastes$associations=='achilles-deidamia' & contrastes$groupes=='Pérou'), 
                            which(contrastes$associations=='helenor-deidamia' & contrastes$groupes=='Pérou'), 
                            which(contrastes$associations=='achilles-achilles' & contrastes$groupes=='Pérou'), 
                            which(contrastes$associations=='helenor-helenor' & contrastes$groupes=='Pérou'), 
                            which(contrastes$associations=='deidamia-deidamia' & contrastes$groupes=='Pérou'), 
                            which(contrastes$associations=='achilles-helenor' & contrastes$groupes=='Guyane'), 
                            which(contrastes$associations=='achilles-deidamia' & contrastes$groupes=='Guyane'),
                            which(contrastes$associations=='helenor-deidamia' & contrastes$groupes=='Guyane'),
                            which(contrastes$associations=='achilles-achilles' & contrastes$groupes=='Guyane'),
                            which(contrastes$associations=='helenor-helenor' & contrastes$groupes=='Guyane'),
                            which(contrastes$associations=='deidamia-deidamia' & contrastes$groupes=='Guyane'),
                            which(contrastes$associations=='achilles-helenor' & contrastes$groupes=='Interlocalités'),
                            which(contrastes$associations=='achilles-deidamia' & contrastes$groupes=='Interlocalités'),
                            which(contrastes$associations=='helenor-deidamia' & contrastes$groupes=='Interlocalités'),
                            which(contrastes$associations=='achilles-achilles' & contrastes$groupes=='Interlocalités'),
                            which(contrastes$associations=='helenor-helenor' & contrastes$groupes=='Interlocalités'),
                            which(contrastes$associations=='deidamia-deidamia' & contrastes$groupes=='Interlocalités')
                            ), ])
  contrastes[ , i] <- apply(contrastes[ , i], 2, function(x) as.numeric(as.character(x)))
  a = ggplot(contrastes) + 
    stat_sum(aes(x=associations, y=dS.mean, ymin=dS.lwr, ymax=dS.upr, col=Observateur), position=position_dodge(0.6), geom = "pointrange", shape=20, size=1) + 
    geom_abline(intercept=1,slope=0,linetype='dashed')+
    ylab('Contraste chromatique (en JND)')+
    xlab('')+
    ggtitle('')+
    facet_grid(~factor(groupes, levels = c('Pérou', 'Guyane', 'Interlocalités')), scales='free_x')+
    theme_linedraw()+
    theme(legend.position = 'none', axis.text = element_text(size=10), axis.title.y = element_text(size = 14))
  b = ggplot(contrastes) + 
    stat_sum(aes(x=associations, y=dL.mean, ymin=dL.lwr, ymax=dL.upr, col=Observateur), position=position_dodge(0.6), geom = "pointrange", shape=20, size=1) + 
    geom_abline(intercept=1,slope=0,linetype='dashed')+
    ylab('Contraste achromatique (en JND)')+
    xlab('')+
    ggtitle('')+
    facet_grid(~factor(groupes, levels = c('Pérou', 'Guyane', 'Interlocalités')), scales='free_x')+
    theme_linedraw()+
    theme(legend.position = 'none', axis.text = element_text(size=10), axis.title.y = element_text(size = 14))
  grid_arrange_shared_legend(a,b, ncol=1, nrow=2, position='bottom')
}

dL.mean = dL.upr = dL.upr = dL.lwr = rep('NA',15)
entrespetpopH = cbind(entrespetpopH, dL.mean, dL.lwr, dL.upr)
graphe(w=entrespetpop, wname='oiseauUV', x=entrespetpopV, xname='oiseauV', y=H_erato_dist_e, yname='H_erato', z=E_atala_dist_e, zname='E_atala', h=entrespetpopH, hname = 'Humains')

#### On récup les femelles ####

setwd('~/Desktop/Femelles')

# Localisation des fichiers 
specfiles = list.files(path='.', pattern='*\\.RFL8',full.names = TRUE, recursive = TRUE)
# Récupération des infos 
infos_fem_sym = do.call(rbind, lapply(specfiles, recup_infos))
infos_fem_sym = infos_fem_sym[which(infos_fem_sym$couleur=="bleu.RFL8"),]
# Importation des spectres de 300 à 700 nm
spectra_fem = getspec("~/Desktop/Femelles", ext = "RFL8", lim = c(300,700), subdir = F)
# Que les bleus 
femel_sympat = subset(spectra_fem, 'bleu')

# Smooth de tous les spectres avec le même coefficient 
spectra_fem_sm = procspec(femel_sympat, opt = 'smooth', fixneg = 'zero', span = 0.2)


Perou_M = subset(symp_males, 'CL')

#spectres
Perou_all = merge(Perou_M, spectra_fem_sm)
achilles_perou= subset(Perou_all, 'achilles')
helenor_perou= subset(Perou_all, 'helenor')
deidamia_perou= subset(Perou_all, 'deidamia')
#infos
infos_perou = rbind(symp_bleu[1:31,1:5], infos_fem_sym)
sp_sex = paste(infos_perou$espece, infos_perou$sex)
achilles = infos_perou[which(infos_perou$espece=='achilles'),]
helenor = infos_perou[which(infos_perou$espece=='helenor'),]
deidamia = infos_perou[which(infos_perou$espece=='deidamia'),]

aggplot(Perou_all, by=sp_sex, legend=T)
par(mfrow=c(2,2))
aggplot(achilles_perou, by=achilles$sex, main='Achilles', legend=T, ylim=c(0,100))
aggplot(helenor_perou, by=helenor$sex, main='Helenor', legend=T, ylim=c(0,100))
aggplot(deidamia_perou, by=deidamia$sex, main='Deidamia', legend=T, ylim=c(0,100))


#ACP 
# Transposition wavelength en variables de l'ACP
preou = t(Perou_all)
# Wavelength en noms de variables 
colnames(preou) = preou[1,]
# Suppression de la colonne WL
preou=preou[-1,]
proeu = prcomp(preou, scale. = T)
summary(proeu)
fviz_pca_ind(proeu, asp=1, labels = F, pointsize=0) + geom_point(aes(color=infos_perou$espece, shape=infos_perou$sex), size=2) + scale_shape_manual(values=c(16,21))
fviz_pca_ind(ACP, axes = c(1,3), asp=1, labels = F, pointsize=0) + geom_point(aes(color=symp_bleu$Loc, shape=symp_bleu$espece), size=2) + scale_shape_manual(values=c(15, 16, 17))


# MANOVA
summary(manova(proeu$x[,1:7]~infos_perou$sex))
summary(manova(proeu$x[,1:7]~infos_perou$espece))
summary(manova(proeu$x[,1:7]~infos_perou$sex*infos_perou$espece))


#### Tests vision models ####

# Modèle oiseaux UV-sensibles
visualmodUV = vismodel(Perou_all, visual = 'avg.uv', achromatic='bt.dc', illum = 'forestshade')
espacecolUV = colspace(visualmodUV, space='tcs')
# Chromatic contrast UV-sensitive birds
entrespetpop = bootcoldist(visualmodUV, n=c(0.37,0.7,0.99,1), weber=0.05, weber.achro=0.05, by=sp_sex)
entrespetpop
graphe_achro(entrespetpop, 'Oiseaux UV')
graphe_chro(entrespetpop, 'Oiseaux UV')

# Modèle oiseaux V-sensibles
visualmodV = vismodel(Perou_all, visual = 'avg.v', achromatic='ch.dc', illum = 'forestshade') #illum = 'forestshade'
espacecolV = colspace(visualmodV, space='tcs')
# Chromatic contrast V-sensitive birds
entrespetpopV = bootcoldist(visualmodV, n=c(0.25,0.5,1,1),weber=0.05,weber.achro=0.5, by=sp_sex)
entrespetpopV
graphe_achro(entrespetpopV, 'Oiseaux V')
graphe_chro(entrespetpopV, 'Oiseaux V')

#### Heliconius visual system ####
## D'après Finkbeiner et al., 2017 ##
H_erato_sens = sensmodel(c(355, 390, 470, 555))
H_erato_vis = vismodel(Perou_all, visual = H_erato_sens, relative = F, achromatic ='all', illum = 'forestshade')
H_erato_colspace = colspace(H_erato_vis, space='tcs')
H_erato_dist_e = bootcoldist(H_erato_vis, by=sp_sex, n=c(0.09, 0.07, 0.17, 1), weber=0.05, weber.achro=0.05) 
H_erato_dist_e
graphe_achro(H_erato_dist_e, 'erato')
graphe_chro(H_erato_dist_e, 'erato')

#### E. atala (Lycaenidae) visual system ####
## D'après Liénard et al., 2021 ##
E_atala_sens = sensmodel(c(362, 441, 495, 564))
E_atala_vis = vismodel(Perou_all, visual = E_atala_sens, relative = F, achromatic = 'all', illum = 'forestshade')
E_atala_colspace = colspace(E_atala_vis, space='tcs')
E_atala_dist_e = bootcoldist(E_atala_vis, by=sp_sex, n=c(0.06, 0.16, 0.03, 0.75), weber=0.05, weber.achro=0.05) #Finkbeiner et al., 2017
E_atala_dist_e
graphe_achro(E_atala_dist_e, 'atala')
graphe_chro(E_atala_dist_e, 'atala')

#### Chromatic contrast Homo sapiens ####
visualmodH = vismodel(Perou_all, visual = 'cie2')
espacecolH = colspace(visualmodH, space='tri')
entrespetpopH = bootcoldist(espacecolH, by=sp_sex)
entrespetpopH
graphe_chro(entrespetpopH, 'humains')

#### Fonction pour graphe tout ensemble ###
graphe = function(w, wname, x, xname, y, yname, z, zname, h, hname){
  lignes = c(1,10,15)
  associations = c('achilles', 'deidamia', 'helenor')
  i = c(1:6)
  tableau_w = w[lignes,]
  tableau_w = as.data.frame(cbind(tableau_w, associations, Observateur = as.character(rep(wname,length(lignes)))))
  tableau_x = x[lignes,]
  tableau_x = as.data.frame(cbind(tableau_x, associations, Observateur = as.character(rep(xname,length(lignes)))))
  tableau_y = y[lignes,]
  tableau_y = as.data.frame(cbind(tableau_y, associations, Observateur = as.character(rep(yname,length(lignes)))))
  tableau_z = z[lignes,]
  tableau_z = as.data.frame(cbind(tableau_z, associations, Observateur = as.character(rep(zname,length(lignes)))))
  tableau_h = h[lignes,]
  tableau_h = as.data.frame(cbind(tableau_h, associations, Observateur = as.character(rep(hname,length(lignes)))))
  contrastes = as.data.frame(rbind(tableau_w, tableau_x, tableau_y, tableau_z, tableau_h))
  contrastes[ , i] <- apply(contrastes[ , i], 2, function(x) as.numeric(as.character(x)))
  a = ggplot(contrastes) + 
    stat_sum(aes(x=associations, y=dS.mean, ymin=dS.lwr, ymax=dS.upr, col=Observateur), position=position_dodge(0.6), geom = "pointrange", shape=20, size=1) + 
    geom_abline(intercept=1,slope=0,linetype='dashed')+
    ylab('Contraste chromatique (en JND)')+
    xlab('')+
    ggtitle('')+
    theme_linedraw()+
    theme(legend.position = 'none', axis.text = element_text(size=12), axis.title.y = element_text(size = 14))
  b = ggplot(contrastes) + 
    stat_sum(aes(x=associations, y=dL.mean, ymin=dL.lwr, ymax=dL.upr, col=Observateur), position=position_dodge(0.6), geom = "pointrange", shape=20, size=1) + 
    geom_abline(intercept=1,slope=0,linetype='dashed')+
    ylab('Contraste achromatique (en JND)')+
    xlab('')+
    ggtitle('')+
    theme_linedraw()+
    theme(legend.position = 'none', axis.text = element_text(size=12), axis.title.y = element_text(size = 14))
  grid_arrange_shared_legend(a,b, ncol=1, nrow=2, position='bottom') 
}

dL.mean = dL.upr = dL.upr = dL.lwr = rep('NA',15)
entrespetpopH = cbind(entrespetpopH, dL.mean, dL.lwr, dL.upr)
graphe(w=entrespetpop, wname='oiseauUV', x=entrespetpopV, xname='oiseauV', y=H_erato_dist_e, yname='H_erato', z=E_atala_dist_e, zname='E_atala', h=entrespetpopH, hname = 'Humains')







