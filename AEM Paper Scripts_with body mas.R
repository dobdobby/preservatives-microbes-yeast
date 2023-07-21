rm(list=ls())

library(MASS)
library(effectsize)
library(car)
library(emmeans)
library(ggplot2)

#setwd("/Volumes/thirdeyesqueegee/Dobslab Dropbox/SanninoDave/TAG yeast and preservatives")

setwd("~/Dobslab Dropbox/SanninoDave/TAG yeast and preservatives")

######################
##	experiment 1	##
##	ax vs cv TAG	##
######################

read.csv("9-17-21 combined.csv",stringsAsFactors = TRUE)->Exp1
  
  #REMOVE 2 MOST EXTREME VALUES
nrow(Exp1)
Exp1 <- subset(Exp1, TAG_mg %in% range(TAG_mg)==F)
nrow(Exp1)

	#normalise TAG to the mean of axenic, no preservatives, yeast S4707
Exp1_normalizer <- mean(subset(Exp1, Bacteria=="Ax" & Yeast=="S4707" & Preservative=="N")$TAG_mg)
Exp1$normalizer <- Exp1_normalizer
Exp1$TAG_mg <- with(Exp1, TAG_mg / normalizer)

	#normalise mg to the mean of axenic, no preservatives, yeast S4707
Exp1$mgFly <- with(Exp1, Weight / number.of.flies)
Exp1_normalizer_mg <- mean(subset(Exp1, Bacteria=="Ax" & Yeast=="S4707" & Preservative=="N")$mgFly)
Exp1$normalizer_mg <- Exp1_normalizer_mg
Exp1$mgFly <- with(Exp1, mgFly / normalizer_mg)

	#factor recoding
levels(Exp1$Yeast) <- c("Yeast A", "Yeast B")
levels(Exp1$Preservative) <- c("Set 1", "Set 2", "None")
Exp1$Preservative <- factor(Exp1$Preservative, levels=c("None", "Set 1", "Set 2"))

	#plot
Exp1plotv = ggplot(Exp1,aes(x=Bacteria,y=TAG_mg, fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha = 0.75) +
  scale_fill_manual(values=c("black","dodgerblue")) +
  labs (title = "Conventional/axenic TAG", x = "", y = "Relative TAG") +
  facet_grid(Yeast~Preservative) +
  theme_bw() +
  theme(legend.position="none")
Exp1plotv

	#model
m_interaction1 <- lm(
  TAG_mg ~ Bacteria*Yeast*Preservative,
  data = Exp1,
  contrasts = list(
    Yeast = "contr.sum",
    Preservative = "contr.sum",
    Bacteria = "contr.sum"
  )
)
	
	#model testing
summary(m_interaction1)
Anova(m_interaction1, type=3)
emmeans_interaction1 <- emmeans(m_interaction1, ~ Bacteria * Yeast * Preservative)
joint_tests(emmeans_interaction1, by=c("Yeast", "Preservative"))
pairs(emmeans_interaction1, by=c("Yeast", "Preservative"))
colsList_exp1 <- list(col=rep(c("black","dodgerblue"), 6))
emmip_interaction1 <- emmip(emmeans_interaction1, 
	~ Bacteria | Yeast | Preservative, CIs=T, 
	dotarg=colsList_exp1,
	CIarg=colsList_exp1)
emmip_interaction1 + 
	theme_bw() +
	labs(title="TAG", x="", y="Relative TAG")
	
	#effect sizes
m_interaction1_ES <- eta_squared(m_interaction1, alternative="t")
m_interaction1_ESplot <- plot(m_interaction1_ES)
m_interaction1_ESplot
m_interaction1_ESplot + 
	theme_bw() +
	labs(title="TAG effect size")

	##body weight##
Exp1plotv_mg = ggplot(Exp1,aes(x=Bacteria,y=mgFly, fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("black","dodgerblue")) +
  labs (title = "Body weight", x = "treatment", y = "mg/fly") +
  facet_grid(Yeast~Preservative) +
  theme_bw()
Exp1plotv_mg

	#model
m_interaction1_mg <- lm(
  mgFly ~ Bacteria*Yeast*Preservative,
  data = Exp1,
  contrasts = list(
    Yeast = "contr.sum",
    Preservative = "contr.sum",
    Bacteria = "contr.sum"
  )
)
	
	#model testing
summary(m_interaction1_mg)
Anova(m_interaction1_mg, type=3)
emmeans_interaction1_mg <- emmeans(m_interaction1_mg, ~ Bacteria * Yeast * Preservative)
joint_tests(emmeans_interaction1_mg, by=c("Yeast", "Preservative"))
pairs(emmeans_interaction1_mg, by=c("Yeast", "Preservative"))
emmip_interaction1_mg <- emmip(emmeans_interaction1_mg, 
	~ Bacteria | Yeast | Preservative, CIs=T, 
	dotarg=colsList_exp1,
	CIarg=colsList_exp1)
emmip_interaction1_mg + 
	theme_bw() +
	labs(title="Body weight", x="", y="Weight mg/fly")

	#effect sizes
m_interaction1_mg_ES <- eta_squared(m_interaction1_mg, alternative="t")
m_interaction1_mg_ESplot <- plot(m_interaction1_mg_ES)
m_interaction1_mg_ESplot
m_interaction1_mg_ESplot + 
	theme_bw() +
	labs(title="Body weight effect size")

######################
##	experiment 2		##
##	gn vs ax TAG		##
######################

read.csv("TAG assay yeast and preservatives 11-9-21 reexamined.csv",stringsAsFactors = TRUE)->Exp2
Exp2$Bacteria <- factor(Exp2$Bacteria, levels=c("Ax","Lb","Ap"))

	#normalise TAG to the mean of axenic, no preservatives, yeast S4707
Exp2_normalizer <- mean(subset(Exp2, Bacteria=="Ax" & Yeast=="S4707" & Preservatives=="N")$TAG_mg)
Exp2$normalizer <- Exp2_normalizer
Exp2$TAG_mg <- with(Exp2, TAG_mg / normalizer)

	#normalise mgFly to the mean of axenic, no preservatives, yeast S4707
Exp2$mgFly <- with(Exp2, weight / number.of.flies)
Exp2_normalizer_mg <- mean(subset(Exp2, Bacteria=="Ax" & Yeast=="S4707" & Preservatives=="N")$mgFly)
Exp2$normalizer_mg <- Exp2_normalizer_mg
Exp2$mgFly <- with(Exp2, mgFly / normalizer_mg)

	#recoding factors
levels(Exp2$Yeast) <- c("Yeast B", "Yeast A")
Exp2$Yeast <- factor(Exp2$Yeast, levels=c("Yeast A", "Yeast B"))
levels(Exp2$Preservatives) <- c("Set 1", "Set 2", "none")
Exp2$Preservatives <- factor(Exp2$Preservatives, levels=c("none", "Set 1", "Set 2"))

	#plot
Exp2plotv = ggplot(Exp2,aes(x=Bacteria,y=TAG_mg, fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("black","purple","red")) +
  labs (title = "TAG", x = "treatment", y = "Relative TAG") +
  facet_grid(Yeast ~ Preservatives) +
  theme_bw() +
  theme(legend.position="none")
Exp2plotv

	#model
m_interaction2 <- lm(
  TAG_mg ~ Bacteria*Yeast*Preservatives,
  data = Exp2,
  contrasts = list(
    Yeast = "contr.sum",
    Preservatives = "contr.sum",
    Bacteria = "contr.sum"
  )
)

	#model testing
summary(m_interaction2)
Anova(m_interaction2, type=3)
emmeans_interaction2 <- emmeans(m_interaction2, ~ Bacteria * Yeast * Preservatives)	

joint_tests(emmeans_interaction2, by=c("Yeast", "Preservatives"))
	#there is a bigger effect of bacterial variation on set 2
pairs(emmeans_interaction2, by=c("Yeast", "Preservatives"))
	#Ap reduces TAG more on set 2, relative to the others
joint_tests(emmeans_interaction2, by=c("Bacteria", "Yeast"))	
	#the overall effect of preservative variation is still significant when Ap present, but diminished (lower F ratios)
pairs(emmeans_interaction2, by=c("Bacteria", "Yeast"))
	#Ap flips the sign of adding set 2 (relative to none) - we go from a sig +ve to sig -ve.

colsList_exp2 <- list(col=rep(c("black","purple","red"), 6))
emmip_interaction2 <- emmip(emmeans_interaction2, 
	~ Bacteria | Yeast | Preservatives, CIs=T, 
	dotarg=colsList_exp2,
	CIarg=colsList_exp2)
emmip_interaction2 + 
	theme_bw() +
	labs(title="TAG", x="", y="Relative TAG (Est. marginal mean ± CI95)")

	#effect sizes
m_interaction2_ES <- eta_squared(m_interaction2, alternative="t")
m_interaction2_ESplot <- plot(m_interaction2_ES)
m_interaction2_ESplot
m_interaction2_ESplot + 
	theme_bw() +
	labs(title="TAG effect size")

	##body mass##
	
	#plot
Exp2plotv_mg = ggplot(Exp2,aes(x=Bacteria,y=mgFly, fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("black","purple","red")) +
  labs (title = "mg/fly", x = "treatment", y = "mg/fly") +
  facet_grid(Yeast ~ Preservatives) +
  theme_bw()
Exp2plotv_mg

	#model
m_interaction2_mg <- lm(
  mgFly ~ Bacteria*Yeast*Preservatives,
  data = Exp2,
  contrasts = list(
    Yeast = "contr.sum",
    Preservatives = "contr.sum",
    Bacteria = "contr.sum"
  )
)

	#model testing
summary(m_interaction2_mg)
Anova(m_interaction2_mg, type=3)
emmeans_interaction2_mg <- emmeans(m_interaction2_mg, ~ Bacteria * Yeast * Preservatives)	
joint_tests(emmeans_interaction2_mg, by=c("Yeast", "Preservatives"))
pairs(emmeans_interaction2_mg, by=c("Yeast", "Preservatives"))	
joint_tests(emmeans_interaction2_mg, by="Bacteria")	
emmip_interaction2_mg <- emmip(emmeans_interaction2_mg, 
	~ Bacteria | Yeast | Preservatives, CIs=T, 
	dotarg=colsList_exp2,
	CIarg=colsList_exp2)
emmip_interaction2_mg + 
	theme_bw() +
	labs(title="Body weight", x="", y="Weight mg/fly")
	
	#effect sizes
m_interaction2_mg_ES <- eta_squared(m_interaction2_mg, alternative="t")
m_interaction2_mg_ESplot <- plot(m_interaction2_mg_ES)
m_interaction2_mg_ESplot
m_interaction2_mg_ESplot + 
	theme_bw() +
	labs(title="Body weight effect size")

############################################################
## compare axenics in experiment 1 and 2  - consistent?	##
############################################################
#all(colnames(Exp1) %in% colnames(Exp2))
#colnames(Exp1) <- tolower(colnames(Exp1))
#colnames(Exp2) <- tolower(colnames(Exp2))
#colnames(Exp2)[6] <- "preservative"
#colnames(Exp2)[14] <- "volume.homogenized"
#all(colnames(Exp1) %in% colnames(Exp2))
#Exp2 <- Exp2[,match(colnames(Exp2), colnames(Exp1))]

#axenics <- rbind(data.frame(subset(Exp1, bacteria=="Ax"), experi = 1),
#data.frame(subset(Exp2, bacteria=="Ax"), experi = 2))
#axenics$experi <- factor(axenics$experi)
#levels(axenics$yeast)[3] <- levels(axenics$yeast)[2]

#AxenicsPlot <- ggplot(axenics,aes(x=experi,y=tag_mg, fill=bacteria)) +
#  geom_violin() +
#  geom_point(col="grey", alpha=0.75) +
#  scale_fill_manual(values=c("black","dodgerblue")) +
#  labs (title = "TAG", x = "treatment", y = "TAG µg/mg") +
#  facet_grid(preservative~yeast) +
#  theme_bw()
#AxenicsPlot

#m_interaction_axenics <- lm(
#  tag_mg ~ experi*yeast*preservative,
#  data = axenics,
#  contrasts = list(
#    yeast = "contr.sum",
#    preservative = "contr.sum",
#    experi = "contr.sum"
#  )
#)

#Anova(m_interaction_axenics, type="III")
#joint_tests(m_interaction_axenics, by="preservative")

######################
##	experiment 3		##
##	pilot CFU quant	##
##	on yeasts C & D	##
######################

read.csv("cfu_fly_log 5-2-23 2.csv",stringsAsFactors = TRUE)->cfufly_pilot_all
cfufly_pilot_all$Bacteria <- factor(cfufly_pilot_all$Bacteria, levels=c("Ax","Lb","Ap"))

str(cfufly_pilot_all)

	#subset without axenics (not informative for modelling)
cfufly_pilot <- droplevels(subset(cfufly_pilot, Bacteria!="Ax"))

	#plot on log scale
cfufly_pilotplotViolin_logScale <- ggplot(cfufly_pilot_all,aes(x=Bacteria, y=CFU_fly, fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("red","black","purple")) +
  labs (title = "CFU/fly", x = "treatment", y = "CFU/fly") +
  theme_bw() +
  facet_grid(Yeast ~ Preservatives) +
  scale_y_continuous(trans="log10")
cfufly_pilotplotViolin_logScale

	#plot log values
cfufly_pilotplotViolin_log <- ggplot(cfufly_pilot_all,aes(x=Bacteria,y=log10(CFU_fly+1), fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("black","purple","red")) +
  labs (title = "CFU/fly", x = "treatment", y = "log10(CFU/fly +1)") +
  theme_bw() +
  facet_grid(Yeast ~ Preservatives)
cfufly_pilotplotViolin_log

#test using a negative binomial model
#using just Ap and Lb data
CFU.nb_pilot <- glm.nb(as.numeric(CFU_fly)  ~ Yeast * Preservatives * Bacteria, data=cfufly_pilot,
                 contrasts = list(
                   Yeast = "contr.sum",
                   Preservatives = "contr.sum",
                   Bacteria = "contr.sum"
                 ), control=glm.control(maxit=1000)
)

Anova(CFU.nb_pilot, type=3)
joint_tests(CFU.nb_pilot)
joint_tests(CFU.nb_pilot, by="Bacteria")

######################
##	experiment 4		##
##	CFU quant + TAG	##
##	on yeasts C&D&E	##
######################

	#TAG data
read.csv("TAG assay yeast and preservatives total data 6-21-23.csv",stringsAsFactors = TRUE) -> Exp4	#SIGMA
#read.csv("TAG assay yeast and preservatives homebrew 6-15-23.csv",stringsAsFactors = TRUE) -> Exp4	#homebrew

str(Exp4)
	
	#reorder factors
Exp4$Bacteria <- factor(Exp4$Bacteria, levels=c("Ax", "Lb", "Ap"))
levels(Exp4$Yeast) <- c("Yeast C", "Yeast D", "Yeast E")
levels(Exp4$Preservatives) <- c("Set 1", "Set 2", "None")
Exp4$Preservatives <- factor(Exp4$Preservatives, levels=c("None", "Set 1", "Set 2"))

	#normalise TAG to the mean of axenic, no preservatives, yeast C
Exp4_normalizer <- mean(subset(Exp4, Bacteria=="Ax" & Yeast=="Yeast D" & Preservatives=="None")$TAG_mg)
Exp4$normalizer <- Exp4_normalizer
Exp4$TAG_mg <- with(Exp4, TAG_mg / normalizer)

	#normalise TAG to the mean of axenic, no preservatives, yeast C
Exp4$mgFly <- with(Exp4, Weight / number.of.flies)
Exp4_normalizer_mg <- mean(subset(Exp4, Bacteria=="Ax" & Yeast=="Yeast D" & Preservatives=="None")$mgFly)
Exp4$normalizer_mg <- Exp4_normalizer_mg
Exp4$mgFly <- with(Exp4, mgFly / normalizer_mg)

	#plot
Exp4plotv = ggplot(Exp4,aes(x=Bacteria,y=TAG_mg, fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("black","purple","red")) +
  labs (title = "TAG", x = "treatment", y = "Relative TAG") +
  facet_grid(Yeast~Preservatives) +
  theme_bw() +
  ylim(c(0.5, 1.5))
Exp4plotv

	#model	
m_interaction4 <- lm(
  TAG_mg ~ Bacteria*Yeast*Preservatives,
  data = Exp4,
  contrasts = list(
    Yeast = "contr.sum",
    Preservatives = "contr.sum",
    Bacteria = "contr.sum"
  )
)

	#test
summary(m_interaction4)
Anova(m_interaction4, type=3)
emmeans_interaction4 <- emmeans(m_interaction4, ~ Bacteria * Yeast * Preservatives)	
joint_tests(emmeans_interaction4, by=c("Yeast", "Preservatives"))
pairs(emmeans_interaction4, by=c("Yeast", "Preservatives"))
joint_tests(emmeans_interaction4, by="Yeast")
joint_tests(emmeans_interaction4, by="Bacteria")	
pairs(emmeans_interaction4, by=c("Yeast","Bacteria"))
colsList_exp4 <- list(col=rep(c("black","purple","red"), 9))
emmip_interaction4 <- emmip(emmeans_interaction4, 
	~ Bacteria | Yeast | Preservatives, CIs=T, 
	dotarg=colsList_exp4,
	CIarg=colsList_exp4)
emmip_interaction4 + 
	theme_bw() +
	labs(title="TAG", x="", y="Relative TAG")

	#effect size
m_interaction4_ES <- eta_squared(m_interaction4, alternative="t")
m_interaction4_ESplot <- plot(m_interaction4_ES)
m_interaction4_ESplot
m_interaction4_ESplot + 
	theme_bw()+
	labs(title="TAG effect size")

	##body mass##

	#plot
Exp4plotv_mg = ggplot(Exp4,aes(x=Bacteria,y=mgFly, fill=Bacteria)) +
  geom_boxplot() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("black","purple","red")) +
  labs (title = "mg/fly", x = "treatment", y = "mg/fly") +
  facet_grid(Yeast ~ Preservatives) +
  theme_bw()
Exp4plotv_mg

	#model
m_interaction4_mg <- lm(
  mgFly ~ Bacteria*Yeast*Preservatives,
  data = Exp4,
  contrasts = list(
    Yeast = "contr.sum",
    Preservatives = "contr.sum",
    Bacteria = "contr.sum"
  )
)

	#test
summary(m_interaction4_mg)
Anova(m_interaction4_mg, type=3)
emmeans_interaction4_mg <- emmeans(m_interaction4_mg, ~ Bacteria * Yeast * Preservatives)	
joint_tests(emmeans_interaction4_mg, by=c("Yeast", "Preservatives"))
pairs(emmeans_interaction4_mg, by=c("Yeast", "Preservatives"))
joint_tests(emmeans_interaction4_mg, by="Yeast")
joint_tests(emmeans_interaction4_mg, by="Bacteria")	
pairs(emmeans_interaction4_mg, by=c("Yeast","Bacteria"))
emmip_interaction4_mg <- emmip(emmeans_interaction4_mg, 
	~ Bacteria | Yeast | Preservatives, CIs=T, 
	dotarg=colsList,
	CIarg=colsList)
emmip_interaction4_mg + 
	theme_bw() +
	labs(title="Body weight", x="", y="Weight mg/fly")
	
	#effect sizes
m_interaction4_mg_ES <- eta_squared(m_interaction4_mg, alternative="t")
m_interaction4_mg_ESplot <- plot(m_interaction4_mg_ES)
m_interaction4_mg_ESplot
m_interaction4_mg_ESplot + 
	theme_bw()+
	labs(title="Body weight effect size")

	#CFU data
CFUdata <- read.csv("TAG assay yeast and preservatives CFUs.csv", stringsAsFactors=T)
CFUdata$Bacteria <- factor(CFUdata$Bacteria, levels=c("Ax", "Lb", "Ap"))

	#relevel
levels(CFUdata$Yeast) <- c("Yeast C", "Yeast D", "Yeast E")
levels(CFUdata$Preservatives) <- c("Set 1", "Set 2", "None")
CFUdata$Preservatives <- factor(CFUdata$Preservatives, levels=c("None", "Set 1", "Set 2"))

	#plot CFUs
CFUplot1 <- ggplot(CFUdata,aes(x=Bacteria,y=CFU_fly, fill=Bacteria)) +
  geom_violin() +
  geom_jitter(col="grey", alpha=0.75, width=0.1, height=0) +
  scale_fill_manual(values=c("black","purple","red")) +
  labs (title = "CFU", x = "") +
  facet_grid(Yeast~Preservatives) +
  theme_bw() +
  ylab("CFU (log10)") +
  theme(legend.position="none")
CFUplot1

	#plot CFUs without axenis
CFUdata_noAx <- droplevels(subset(CFUdata, Bacteria!="Ax"))
CFUplot2 = ggplot(CFUdata_noAx,aes(x=Bacteria,y=CFU_fly, fill=Bacteria)) +
  geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  scale_fill_manual(values=c("purple","red")) +
  labs (title = "CFU", x = "treatment", y = "TAG/mg") +
  facet_grid(Yeast~Preservatives) +
  theme_bw() +
  ylab("CFU / log (log10)")
CFUplot2

	#statistical effect of yeast * preservatives on CFU?
CFUmod <- lm(
  CFU_fly ~ Bacteria*Yeast*Preservatives,
  data = CFUdata_noAx,
  contrasts = list(
    Yeast = "contr.sum",
    Preservatives = "contr.sum",
    Bacteria = "contr.sum"
  )
)

	#test
Anova(CFUmod, type=3)
joint_tests(CFUmod, by="Bacteria")
emmeans_CFU <- emmeans(CFUmod, ~ Bacteria * Yeast * Preservatives)	
pairs(emmeans_CFU, by=c("Yeast", "Preservatives"))
pairs(emmeans_CFU, by=c("Bacteria", "Yeast"))
pairs(emmeans_CFU, by=c("Bacteria", "Preservatives"))

	#calculate mean CFU and mean TAG per condition
meanCFU <- aggregate(CFU_fly ~ Preservatives * Yeast * Bacteria, CFUdata, mean)
meanTAG <- aggregate(TAG_mg ~ Preservatives * Yeast * Bacteria, Exp4, mean)
means <- data.frame(meanCFU, TAG_mg=meanTAG$TAG_mg)
ggplot(means, aes(x=CFU_fly, y=TAG_mg, group=Bacteria, col=Bacteria, shape=Preservatives))+
  #	geom_violin() +
  geom_point(col="grey", alpha=0.75) +
  geom_smooth(method="lm", se=F) 