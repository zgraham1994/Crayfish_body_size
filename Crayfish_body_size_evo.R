#Test run of PCM for Manuscript on Body Size Evolution in Crayfishes
#9/25/23


#Packages
library(ape) #Analysis of Phylogenetics and Evolution in r, "APE"
library(geiger)
library(phytools)
library(phangorn)
library(nlme)
library(beeswarm)
library(scales)
library(emmeans)
library(phylolm)

#Reset R's graphing parameters
#old.par <- par(mar = c(0, 0, 0, 0))
#par(old.par)
###                            ###
#Load in date, Prune tree, etc####
###                            ###
#read in  data
cray.data<-read.csv("Cray_body_size_max.csv", row.names=1, stringsAsFactors = TRUE)
head(cray.data)
str(cray.data)

#read in tree
cray.tree<-read.tree("stern.2017.ml.new.txt")
cray.tree
summary(cray.tree)  

#Check to see if data and phylogeny match
chk<-name.check(cray.tree,cray.data)
summary(chk)
#chk$data_not_tree
#chk$tree_not_data
#They do not, there are taxa in the tree but not in the data, and taxa present in the data
#but not the tree

#Align the tree and data: first with drop.tip, then with prune
tree.pruned = drop.tip(cray.tree,chk$tree_not_data)
data.pruned = cray.data[tree.pruned$tip.label,]
name.check(tree.pruned,data.pruned) #This should return "OK" then we are good to go

###                         ###
#Figure 1: Circle Phylogeny####
###                         ###

#I need to make an ancestral reconstruction first, and then create a simmap
#If i do not do this, I cannot get the colors mapped to tip labels, so it is just 
#an extra step to get the figure looking nice
body.size<-setNames(data.pruned$max.of.cl.and.total.length,rownames(data.pruned))

#Comparing SYM,ER,and ARD models; 1000 simulations each
burrow = setNames(data.pruned$habitat.category,rownames(data.pruned))

bur.ER = fitDiscrete(tree.pruned, burrow, model = "ER", control = list(niter = 1000), ncores = 4)
bur.ER
#model summary:
#log-likelihood = -305.808102
#AIC = 613.616205
#AICc = 613.625094
#free parameters = 1

bur.ARD = fitDiscrete(tree.pruned, burrow, model = "ARD", control = list(niter = 1000), ncores = 4)
bur.ARD
#model summary:
# log-likelihood = -283.416380
#AIC = 590.832759
#AICc = 591.543465
#free parameters = 12

bur.SYM = fitDiscrete(tree.pruned, burrow, model = "SYM", control = list(niter = 1000), ncores = 4)
bur.SYM
#model summary:
#log-likelihood = -291.855416
#AIC = 595.710832
#AICc = 595.899596
#free parameters = 6

#ARD and SYM have AIC within 2 values (BARELY); so both are likely
#SYM is more parsimonious though, so I am going to use that model

#Then I use make.simmap to simulate the state across the tree which I can use to color nodes based on results
#Using ARD because it is the lowest AIC and most parsimonious
#bur.simARD = make.simmap(tree.pruned, burrow, model = "ARD", nsim = 1000, Q = "empirical")

colorbars2 <-setNames(c(alpha('#bc672c',1),alpha('black',1),alpha('#79c5e7',1),alpha('purple',1)), c("Semi-terrestrial", "Semi-terrestrial/Aquatic","Aquatic","Cave"))
ss2 <-getStates(bur.simARD, "tips")
str(ss2)
barcols2 <- setNames(sapply(ss2, function(x, y) y[which(names(y) == x)], y = colorbars2), names(ss2))

colors2 <- c(c(alpha('#79c5e7',1),alpha('purple',1),alpha('#bc672c',1),alpha('black',1)))
cols2 <- setNames(colors2[1:length(unique(burrow))], sort(unique(burrow)))

test <- summary(bur.simARD, plot = FALSE)
test

#Getting the Posterior Probability of Ancestral Nodes
test$ace
nodelabels(cex=.4)#make node labels visible; root is 454; 584 is root of Astacoidea; 454 is root ot Parstacoidea
#           Aquatic  Cave     Semi-terrestrial    Semi-terrestrial/Aquatic
#454         0.513    0.002    0.398               0.087


#Plot Tree.wbars
#png("Fig1.png",res=300,w=200,h=200,units='mm')

data.pruned$spp <- rownames(data.pruned)
unique(data.pruned$family)
Astacidae <- data.pruned$spp[data.pruned$family == "astacidae"]
Cambaridae <- data.pruned$spp[data.pruned$family == "cambaridae"]
Cambaroididae <- data.pruned$spp[data.pruned$family == "cambaroididae"]
Parastacidae <- data.pruned$spp[data.pruned$family == "parastacidae"]

node1 <- getMRCA(tree.pruned, Astacidae)
node2 <- getMRCA(tree.pruned, Cambaridae)
node3 <- getMRCA(tree.pruned, Cambaroididae)
node4 <- getMRCA(tree.pruned, Parastacidae)

nodes <- c(node1, node2, node3, node4)
nodes

labels <- c("Astacidae", "Cambaridae", "Cambaroididae", "Parastacidae")
#lncol <- c("red", "blue", "green", "yellow")
lncol <- c("gray", "black", "black", "gray")

par(fg = "transparent", col.lab = "white", col = "white")

plotTree.wBars(tree.pruned, body.size, type = "fan", lwd = 1.2, tip.labels = TRUE, ftype = "i", fsize = 0.03,
               border = FALSE, col = barcols2)

nodelabels(pie = test$ace, piecol = cols2, cex = 0.2, border = TRUE)

for(i in 1:length(nodes)) 
  arc.cladelabels(text = "", node = nodes[i], ln.offset = 1.02, lab.offset = 1.1, mark.node = FALSE, lwd = 3,
                  col = lncol[i], orientation = if(labels[i] %in% c("Astacidae", "Cambaridae", "Cambaroididae", "Parastacidae"))
                    "horizontal" else "curved")
par(col.lab = "black", col = "black")
#legend("topleft", legend = labels, lwd = 3, lty = 1, col = lncol, bty = "n")

nodelabels(pie = test$ace, piecol = cols2, cex = 0.2)

#Add legend
cols<-setNames(c(alpha('#bc672c',1),alpha('black',1),alpha('#79c5e7',1),alpha('purple',1)), c("Semi-terrestrial", "Semi-terrestrial/Aquatic","Aquatic","Cave"))

leg<-legend(x="top",legend=names(cols), inset=0.52,
            pch=22,pt.bg=cols,pt.cex=2,
            title="Habitat", cex=.9, bty="n")
#dev.off()

#PGLS Analysis for Habitat and Body Size####
#Built the model
model.1<- phylolm(log(max.of.cl.and.total.length)~habitat.category,phy=tree.pruned,data=data.pruned,boot=1000, btol=30)
summary(model.1)
aov(model.1)

###                                          ###
#Figure 2: Boxplot of Habitat and body size####
###                                          ###
#dev.new(width=10, height=5, unit="in")
data.pruned$habitat.category <- factor(data.pruned$habitat.category , levels=c("Cave","Semi-terrestrial","Semi-terrestrial/Aquatic","Aquatic"))
#png("Fig2.png",res=600,w=250,h=150,units='mm')
beeswarm(log(max.of.cl.and.total.length)~habitat.category, data=data.pruned, spacing=.7, pch=19,cex= 1.5,cex.lab=1.8,cex.axis=1.2, xlab = "Habitat", ylab = "Maximum Carapace Length (log)", col=c(alpha('purple',0.75),alpha('#bc672c',0.75),alpha('black',0.75),alpha('#79c5e7',0.75)))
boxplot(log(max.of.cl.and.total.length)~habitat.category, data=data.pruned,cex.axis=0, add = T,outline=FALSE, names = c("","","",""), col=c(alpha("gray",0.25)))
#dev.off() 


###                              ###
#Temperature/Bergmann's Analysis#####
###                              ###
#Loading in separate date with cave species removed
#read in  data
cray.data.nocave<-read.csv("Cray_body_size_max.nocave.csv", row.names=1, stringsAsFactors = TRUE,strip.white = TRUE)
head(cray.data.nocave)
str(cray.data.nocave)

#read in tree
cray.tree<-read.tree("stern.2017.ml.new.txt")
cray.tree
summary(cray.tree)  

#Check to see if data and phylogeny match
chk<-name.check(cray.tree,cray.data.nocave)
summary(chk)
#chk$data_not_tree
#chk$tree_not_data
#They do not, there are taxa in the tree but not in the data, and taxa present in the data
#but not the tree

#Align the tree and data: first with drop.tip, then with prune
tree.pruned.nocave = drop.tip(cray.tree,chk$tree_not_data)
data.pruned.nocave = cray.data.nocave[tree.pruned.nocave$tip.label,]
name.check(tree.pruned.nocave,data.pruned.nocave) #This should return "OK" then we are good to go
 
#PGLS Analysis for Body size, Habitat*Temperature####
#Built the model
model.2 <- summary(phylolm(log(max.of.cl.and.total.length)~tmean.avg*habitat.category,phy=tree.pruned.nocave,data=data.pruned.nocave,boot=1000, btol=30))
model.2

###                                                                             ###
#Figure 3: Scatter plot of Mean Annual Temperature and Carapace Length Separated####
###                                                                             ###
#Subset data
semi.terrestrial<-data.pruned.nocave[ which(data.pruned.nocave$habitat.category=='Semi-terrestrial'), ]
semi.terrestrial.aquatic <- data.pruned.nocave[ which(data.pruned.nocave$habitat.category=='Semi-terrestrial/Aquatic'), ]
aquatic <- data.pruned.nocave[ which(data.pruned.nocave$habitat.category=='Aquatic'), ]

#Semi.terrestrial
#png("Fig3.png",res=600,w=350,h=150,units='mm')
par(mfrow=c(1,3))
par(mar=c(5,6,4,1)+.1)

colors <- c(alpha("#bc672c",0.75))
shapes = c(19)

plot(semi.terrestrial$tmean.avg,log(semi.terrestrial$max.of.cl.and.total.length), 
     pch=shapes, cex=2,bty='l', xaxt="s", 
     col=colors, las=1, cex.lab= 2, xlim = c(0,28),ylim = c(3,6), 
     xlab="Mean Annual Temperature (°C)", ylab="Maximum Carapace Length (log)")

curve(((3.992215 + 0.672367) + ((0.042534-0.050398)*x)), col=alpha('#bc672c',0.75),lwd=3, from=5, to=22, lty=2, add=T) #Primary
legend("topleft","A)",bty='n',cex= 2)

#Semi.terrestrial/Aquatic
colors <- c(alpha("black",0.75))
shapes = c(19)

plot(semi.terrestrial.aquatic$tmean.avg,log(semi.terrestrial.aquatic$max.of.cl.and.total.length), 
     pch=shapes, cex=2,bty='l', xaxt="s", 
     col=colors, las=1, cex.lab= 2, xlim = c(0,28),ylim = c(3,6), 
     xlab="Mean Annual Temperature (°C)", ylab="Maximum Carapace Length (log)")

curve(((3.992215  + 0.872914) + ((0.042534 -0.051654)*x)), col=alpha('black',0.75),lwd=3,from=2, lty=2, add=T) #Secondary
legend("topleft","B)",bty='n',cex= 2)

#Aquatic
colors <- c(alpha("#79c5e7",0.75))
shapes = c(19)

plot(aquatic$tmean.avg,log(aquatic$max.of.cl.and.total.length), 
     pch=shapes, cex=2,bty='l', xaxt="s", 
     col=colors, las=1, cex.lab= 2, xlim = c(0,28),ylim = c(3,6), 
     xlab="Mean Annual Temperature (°C)", ylab="Maximum Carapace Length (log)")

curve(((3.992215) + (0.042534*x)), col=alpha('#79c5e7',0.75),lwd=3, lty=2, add=T) #Tertiary
legend("topleft","C)",bty='n',cex= 2)
#dev.off() 

###                                                                             ###
#Figure 3 w/gray points: Scatter plot of Mean Annual Temperature and Carapace Length Separated####
###                                                                             ###
png("Fig3.withgray.png",res=600,w=350,h=150,units='mm')
#Semi.terrestrial
par(mfrow=c(1,3))
par(mar=c(5,8,4,1)+.1)
par(mgp=c(4,1,0))

colors <- c(alpha("gray",0.3), alpha("#bc672c",0.8), alpha("gray",0.3))
colors <- colors[as.factor(data.pruned.nocave$habitat.category)]

shapes = c(19,19,19)
shapes <- shapes[as.factor(data.pruned.nocave$habitat.category)]

plot(data.pruned.nocave$tmean.avg,log(data.pruned.nocave$max.of.cl.and.total.length), 
     pch=shapes, cex=2,bty='l', xaxt="s", 
     col=colors, las=1, cex.lab= 2,cex.axis = 1.8, xlim = c(0,28),ylim = c(3,6), 
     xlab="Mean Annual Temperature (°C)", ylab="Maximum Carapace Length (log)")

curve(((3.992215 + 0.672367) + ((0.042534-0.050398)*x)), col=alpha('#bc672c',0.75),lwd=3, from=5, to=22, lty=2, add=T) #Primary
legend("topleft","A)",bty='n',cex= 2)

#Semi.terrestrial/Aquatic
colors <- c(alpha("gray",0.3), alpha("gray",0.3), alpha("black",0.8))
colors <- colors[as.factor(data.pruned.nocave$habitat.category)]

shapes = c(19,19,19)
shapes <- shapes[as.factor(data.pruned.nocave$habitat.category)]

plot(data.pruned.nocave$tmean.avg,log(data.pruned.nocave$max.of.cl.and.total.length), 
     pch=shapes, cex=2,bty='l', xaxt="s", 
     col=colors, las=1, cex.lab= 2,cex.axis = 1.8,  xlim = c(0,28),ylim = c(3,6), 
     xlab="Mean Annual Temperature (°C)", ylab="Maximum Carapace Length (log)")

curve(((3.992215  + 0.872914) + ((0.042534 -0.051654)*x)), col=alpha('black',0.75),lwd=3,from=2, lty=2, add=T) #Secondary
legend("topleft","B)",bty='n',cex= 2)


#Aquatic
#par(mfrow=c(1,3))
par(mar=c(5,6,4,1)+.1)

colors <- c(alpha("#79c5e7",0.8), alpha("gray",0.3), alpha("gray",0.3))
colors <- colors[as.factor(data.pruned.nocave$habitat.category)]

shapes = c(19,19,19)
shapes <- shapes[as.factor(data.pruned.nocave$habitat.category)]

plot(data.pruned.nocave$tmean.avg,log(data.pruned.nocave$max.of.cl.and.total.length), 
     pch=shapes, cex=2,bty='l', xaxt="s", 
     col=colors, las=1, cex.lab= 2,cex.axis = 1.8,  xlim = c(0,28),ylim = c(3,6), 
     xlab="Mean Annual Temperature (°C)", ylab="Maximum Carapace Length (log)")

legend("topleft","C)",bty='n',cex= 2)
curve(((3.992215) + (0.042534*x)), col=alpha('#79c5e7',0.75),lwd=3, lty=2, add=T) #Tertiary
dev.off() 

#Phylogenetic Signal (Blombergs K) ####
#Calculating Blombergs K (Phylogenetic Signal; tendency for related species to resemble one another more than
#expected by chance)
phylosig(tree.pruned, body.size)
#Phylogenetic signal K : 0.048072 

#Test for significance
K_bs<-(phylosig(tree.pruned,body.size,
                test=TRUE,nsim=10000))
K_bs
#Phylogenetic signal K : 0.0473317  
#P-value (based on 10000 randomizations) : 0.3825 
#Not significant, lets plot
plot(K_bs, las=1, cex.axis=0.9)
#This demonstrates that species DO NOT resemble one another more than expected by chance (i.e., ecology, constraint?)

#Fitting Evolutionary Model to find ancestral state of root####
#1.23.24
#Now we can fit a BM model using fitContinous from Gieger
fitBM_bs<-fitContinuous(tree.pruned,body.size)
fitBM_bs
#GEIGER-fitted comparative model of continuous data
#fitted 'BM' model parameters:
#sigsq = 423.667670
#z0 = 107.479107

#Compare BM, EB, and OU
fitEB_bs<-fitContinuous(tree.pruned,body.size,
                        model="EB")
fitEB_bs

#fitOU_bs takes an incredibly long time to run
fitOU_bs<-fitContinuous(tree.pruned,body.size,
                        model="OU")
fitOU_bs

aic_bs<-setNames(c(AIC(fitBM_bs),
                   AIC(fitEB_bs),AIC(fitOU_bs)),
                 c("BM","EB","OU"))
aic_bs
aic.w(aic_bs)
#BM       EB       OU 
#5567.052 5569.061 5008.157 

#Weights
#BM EB OU 
#0  0  1

fitBM_bs<-fitContinuous(tree.pruned,body.size)
fitBM_bs
#GEIGER-fitted comparative model of continuous data
#fitted 'BM' model parameters:
#sigsq = 423.667670
#z0 = 107.479107




###Everything below here is old code####
###                                                                             ###
#Figure 3: Scatter plot of Mean Annual Temperature and Carapace Length Combined###
###                                                                             ###
#png("Fig3.png",res=300,w=200,h=200,units='mm')
colors <- c(alpha("#79c5e7",0.75), alpha("#bc672c",0.75), alpha("black",0.75))
#colors <- c(alpha("black",0.75), alpha("black",0.75), alpha("black",0.75), alpha("black",0.75)) #burrow colors
colors <- colors[as.factor(data.pruned.nocave$habitat.category)]

shapes = c(19,19,19)
shapes <- shapes[as.factor(data.pruned.nocave$habitat.category)]

#png("Cambarus_bergmann.png",res=600,w=300,h=240,units='mm')
par(mar = c(5, 5, 4, 1) + 0.1)

plot(data.pruned.nocave$tmean.avg,log(data.pruned.nocave$max.of.cl.and.total.length), 
     pch=shapes, cex=2,bty='l', 
     col=colors, las=1, cex.lab= 2,
     xlab="Mean Annual Temperature (°C)", ylab="Maximum Carapace Length (log)")

curve((( 4.664582) + (-0.007864*x)), col=alpha('#bc672c',0.75),lwd=3, from=5, to=22, lty=2, add=T) #Primary
curve(((4.664582  + 0.200547) + ((-0.007864 -0.001256)*x)), col=alpha('black',0.75),lwd=3,from=2, lty=2, add=T) #Secondary
curve(((4.664582  -0.672367) + ((-0.007864 + 0.050398 )*x)), col=alpha('#79c5e7',0.75),lwd=3, lty=2, add=T) #Tertiary
#dev.off()


#PGLS Analysis for Burrowing and Temperature###
#Built the model
model.3 <- summary(phylolm(tmean.avg~habitat.category,phy=tree.pruned,data=data.pruned.nocave,boot=100, btol=30))
model.3

###                                                    ###
#Figure X: Boxplot of burrowing and  lattitude midpoint####
###                                                    ###

#dev.new(width=10, height=5, unit="in")
#png("FigXX.png",res=600,w=250,h=150,units='mm')
beeswarm(tmean.avg~habitat.category, data=data.pruned.nocave, spacing=.7, pch=19,cex= 1.5, xlab = "Burrowing Classification", ylab = 'Mean Annual Temperature (°C)', col=c(alpha('#bc672c',0.75),alpha('black',0.75),alpha('#79c5e7',0.75)))
boxplot(tmean.avg~habitat.category, data=data.pruned.nocave, add = T,outline=FALSE, names = c("","",""), col=c(alpha("gray",0.25)))
#dev.off()

#Fitting Evolutionary Model##
#1.23.24
#Now we can fit a BM model using fitContinous from Gieger
fitBM_bs<-fitContinuous(tree.pruned,body_size)
fitBM_bs
#GEIGER-fitted comparative model of continuous data
#fitted 'BM' model parameters:
#sigsq = 423.667670
#z0 = 107.479107


#Compare BM, EB, and OU
fitEB_bs<-fitContinuous(tree.pruned,body_size,
                        model="EB")
fitEB_bs

#fitOU_bs takes an incredibly long time to run
fitOU_bs<-fitContinuous(tree.pruned,body_size,
                        model="OU")
fitOU_bs

aic_bs<-setNames(c(AIC(fitBM_bs),
                   AIC(fitEB_bs),AIC(fitOU_bs)),
                 c("BM","EB","OU"))
aic_bs
aic.w(aic_bs)
#BM       EB       OU 
#5579.445 5581.454 5022.637 

