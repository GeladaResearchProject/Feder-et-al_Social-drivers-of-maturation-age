## LOAD PACKAGES NEEDED
require(plyr);require(dplyr);require(lme4); require(lmerTest); require(MuMIn); require(ggplot2);require(ggplot2);require(ggthemes);require(ggbeeswarm);require(effects);require(coxme);require(reshape2)
theme_set(theme_classic(base_size=18))

## DATA CORRELATIONS

unique.fems<-maturation.months[!duplicated(maturation.months$Individual),]

summary(lm(data=unique.fems, UnitSize_3.5~UnitSize_Conception))
plot(data=unique.fems, UnitSize_3.5~UnitSize_Conception)
summary(lm(data=unique.fems, Rank_3.5~Rank_Conception))

## FIT MODEL WITHOUT TIME TRANSFORMATION
maturation.months$End.Age[maturation.months$Start.Age==maturation.months$End.Age]<-maturation.months$End.Age[maturation.months$Start.Age==maturation.months$End.Age] + 0.001
surv_object<-Surv(time = maturation.months$Start.Age/12, time2=maturation.months$End.Age/12, event=maturation.months$Matured)
fit.coxph.no.time <- coxph(surv_object ~ Takeover.Range + Rank_Conception + scale(UnitSize_Conception) + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)

schoenfeld<-cox.zph(fit.coxph.no.time)
schoenfeld
dev.off()
par(mfrow=c(1,2))

## PLOT SCHOENFELD RESIDUALS
plot(schoenfeld[4],lwd=2, ylab="Beta(t) for father presence", xlab="Age (years)", ylim=c(-4,5))
abline(0,0, col=1,lty=3,lwd=2)
abline(h= fit.coxph.no.time$coef[4], col=2, lwd=2, lty=2)
title(outer=T,adj=0,main="a",cex=1.1,col="black",font=1,line=-2)

plot(schoenfeld[3],lwd=2,ylab="Beta(t) for unit size at conception", xlab="Age (years)",ylim=c(-4,5))
abline(0,0, col=1,lty=3,lwd=2)
abline(h= fit.coxph.no.time$coef[3], col=2, lwd=2, lty=2)
title(outer=T,adj=0.5,main="b",cex=1.1,col="black",font=1,line=-2)

## CONSTRUCT FOUR MODELS
fit.coxph1.poly <- coxph(surv_object ~ Takeover.Range + Rank_Conception + tt(UnitSize_Conception) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph2.poly <- coxph(surv_object ~ Takeover.Range + Rank_Conception + tt(UnitSize_3.5) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph3.poly <- coxph(surv_object ~ Takeover.Range + Rank_3.5 + tt(UnitSize_Conception) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph4.poly <- coxph(surv_object ~ Takeover.Range + Rank_3.5 +  tt(UnitSize_3.5) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph.NULL <- coxph(surv_object ~ 1,  data = maturation.months)

MuMIn::model.sel(fit.coxph1.poly, fit.coxph2.poly, fit.coxph3.poly, fit.coxph4.poly, fit.coxph.NULL)
fit.coxph1.poly
fit.coxph2.poly
fit.coxph3.poly
fit.coxph4.poly

## ASSESS QUADRATIC IMPROVEMENT

model1<-coxph(surv_object ~ Takeover.Range + Rank_Conception + UnitSize_Conception + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
model2<-coxph(surv_object ~ Takeover.Range + Rank_Conception + poly(UnitSize_Conception,2) + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
summary(model1)
summary(model2)
model.sel(model1, model2)

model1<-coxph(surv_object ~ Takeover.Range + Rank_Conception + UnitSize_3.5 + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
model2<-coxph(surv_object ~ Takeover.Range + Rank_Conception + poly(UnitSize_3.5,2) + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
summary(model1)
summary(model2)
model.sel(model1, model2)

## NO IMPROVEMENT

## TRY WITH TIME-VARYING UNIT SIZE

fit.coxph.varying <- coxph(surv_object ~ Takeover.Range + Rank_Conception + tt(Varying_Females) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
summary(fit.coxph.varying)

maturation.months$CurrCat<-"SMALL"
maturation.months$CurrCat[maturation.months$Varying_Females>4]<-"MEDIUM"
maturation.months$CurrCat[maturation.months$Varying_Females>7]<-"LARGE"
maturation.months$CurrCat<-factor(maturation.months$CurrCat, levels=c("SMALL", "MEDIUM", "LARGE"))
fit.group.curr <- survfit(surv_object ~ CurrCat, data = maturation.months)

unit.palette<-c("#993404", "#ec7014","#fec44f")
plot.TV<-ggsurvplot(fit.group.curr, data = maturation.months, xlim=c(4,6.5), censor=F, break.time.by=.5, xlab="Age (years)", ylab="Proportion yet to mature", legend.labs=c("Small", "Medium", "Large"), legend.title="", legend="top", palette = unit.palette)  
plot.TV<-plot.TV$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot.TV

## RESULTS CONGRUENT

## ASSESS UNCERTAINTY IN SURVIVAL MODEL

models1<-list()
models2<-list()
P.rank1<-c()
P.unit1<-c()
P.rank2<-c()
P.unit2<-c()

maturations<-maturation.months[maturation.months$Matured==1,]
mat.ranges$Uncertainty<-NA
for(i in 1:1000) {
  for(j in 1:nrow(mat.ranges)) {
    mat.ranges$Uncertainty[j]<-runif(n=1, min=mat.ranges$DownRange[j], max=mat.ranges$UpRange[j])
  }
  sim.data<-merge(maturation.months, mat.ranges, by="Individual")
  sim.data$Start.Age<-round((sim.data$Start.Age+sim.data$Uncertainty/30.4),1)
  sim.data$End.Age<-round((sim.data$End.Age+sim.data$Uncertainty/30.4),1)
  sim.data$End.Age[sim.data$Start.Age==sim.data$End.Age]<- sim.data$End.Age[sim.data$Start.Age==sim.data$End.Age]+0.01
  surv_object<-Surv(time=sim.data$Start.Age/12, time2=sim.data$End.Age/12, event=sim.data$Matured)
  models1[[i]]<-coxph(surv_object ~ Takeover.Range + Rank_Conception + UnitSize_Conception + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = sim.data)
  P.rank1[i]<-summary(models1[[i]])$coefficients[2,6]
  P.unit1[i]<-summary(models1[[i]])$coefficients[3,6]
  models2[[i]]<-coxph(surv_object ~ Takeover.Range + Rank_Conception + UnitSize_3.5 + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = sim.data)
  P.rank2[i]<-summary(models2[[i]])$coefficients[2,6]
  P.unit2[i]<-summary(models2[[i]])$coefficients[3,6]
}

sum(P.rank1>0.05)/1000
sum(P.unit1>0.05)/1000
sum(P.rank2>0.05)/1000
sum(P.unit2>0.05)/1000

summary(model.avg(models1))
summary(model.avg(models2))

## GENERATE PREDICTED CURVES
## RANK CURVES
library(broom)
library(tidyr)
library(magrittr)
surv_object<-Surv(time = maturation.months$Start.Age/12, time2=maturation.months$End.Age/12, event=maturation.months$Matured)

fit.coxph.no.time <- coxph(surv_object ~ Takeover.Range + Rank_Conception + UnitSize_3.5 + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
newdf <- maturation.months %$%
  expand.grid(UnitSize_3.5=6, Rank_Conception=c(0,1), Takeover.Range=0, Dad="N", Cumulative.Rainfall=mean(maturation.months$Cumulative.Rainfall), Avg.Min=mean(maturation.months$Avg.Min), MajorInjury="N")
newdf$group <- as.character(1:nrow(newdf))

survcurv1 <- survfit(fit.coxph.no.time, newdata = newdf) %>%
  tidy() %>% 
  gather('key', 'value', -time, -n.risk, -n.event, -n.censor) %>% 
  mutate(group = substr(key, nchar(key), nchar(key)),
         key   = substr(key, 1, nchar(key) - 2)) %>% 
  left_join(newdf, 'group') %>% 
  spread(key, value)


## LOW: 63.1 months median
## HIGH: 56.5 months median

ggplot(survcurv1, aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                      col = as.character(Rank_Conception), fill = as.character(Rank_Conception))) +
  geom_line(size = 1)+
  geom_ribbon(alpha = 0.1, col = NA) + xlim(3.5,6)
survcurv1$time<-survcurv1$time*12

## UNIT SIZE CURVES
newdf <- maturation.months %$%
  expand.grid(UnitSize_3.5=c(2,6,10), Rank_Conception=0.5, Takeover.Range=0, Dad="N", Cumulative.Rainfall=mean(maturation.months$Cumulative.Rainfall), Avg.Min=mean(maturation.months$Avg.Min), MajorInjury="N")
newdf$group <- as.character(1:nrow(newdf))

survcurv2 <- survfit(fit.coxph.no.time, newdata = newdf) %>%
  tidy() %>% 
  gather('key', 'value', -time, -n.risk, -n.event, -n.censor) %>% 
  mutate(group = substr(key, nchar(key), nchar(key)),
         key   = substr(key, 1, nchar(key) - 2)) %>% 
  left_join(newdf, 'group') %>% 
  spread(key, value)

survcurv2
ggplot(survcurv2, aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                      col = as.character(UnitSize_3.5), fill = as.character(UnitSize_3.5))) +
  geom_line(size = 1)+
  geom_ribbon(alpha = 0.1, col = NA) + xlim(3.5,6)

survcurv2$time<-survcurv2$time*12

## SMALL: 62.6 months median
## MEDIUM: 60.9 months median
## LARGE: 56.8 months median

## GENERATE SURVIVAL CURVES FOR PLOTS
maturation.months$RankCat_Conception<-factor(maturation.months$RankCat_Conception, levels=c("Low", "High"))
maturation.months$Unit_Cat_Conception<-factor(maturation.months$Unit_Cat_Conception, levels=c("Small", "Medium", "Large"))
maturation.months$RankCat_3.5<-factor(maturation.months$RankCat_3.5, levels=c("Low", "High"))
maturation.months$Unit_Cat_3.5<-factor(maturation.months$Unit_Cat_3.5, levels=c("Small", "Medium", "Large"))
maturation.months$Dad<-factor(maturation.months$Dad, levels=c("N","Y"))

maturation.months$UnitCat1<-"Small"
maturation.months$UnitCat1[maturation.months$UnitSize_Conception>6]<-"Large"
maturation.months$UnitCat1<-factor(maturation.months$UnitCat1, levels=c("Small", "Large"))

maturation.months$UnitCat2<-"Small"
maturation.months$UnitCat2[maturation.months$UnitSize_3.5>7]<-"Large"
maturation.months$UnitCat2<-factor(maturation.months$UnitCat2, levels=c("Small", "Large"))

fit.rank.conception <- survfit(surv_object ~ RankCat_Conception, data = maturation.months)
fit.rank.3.5 <- survfit(surv_object ~ RankCat_3.5, data = maturation.months)
fit.unit.conception <- survfit(surv_object ~ Unit_Cat_Conception, data = maturation.months)
fit.unit.3.5 <- survfit(surv_object ~ Unit_Cat_3.5, data = maturation.months)
fit.fathers <- survfit(surv_object ~ Dad, data = maturation.months)
fit.rank.inj <- survfit(surv_object ~ MajorInjury, data = maturation.months)
fit.1 <- survfit(surv_object ~ UnitCat1, data = maturation.months)
fit.2 <- survfit(surv_object ~ UnitCat2, data = maturation.months)

rank.palette<-c("#440154FF","#21908CFF")
unit.palette<-c("#993404", "#ec7014","#fec44f")
dad.palette<-c("#1F78B4", "#33A02C")
inj.palette<-scale_color_brewer(palette="Dark2")

library(survminer)
plot1a<-ggsurvplot(fit.rank.conception, data = maturation.months, xlim=c(4,6.5), censor=F, break.time.by=.5, xlab="Age (years)", ylab="Proportion yet to mature", legend.labs=c("Low-ranking", "High-ranking"), legend.title="", legend="top", palette = rank.palette)  
plot1a<-plot1a$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1a
plot1b<-ggsurvplot(fit.unit.3.5, data = maturation.months, xlim=c(4,6.5), legend="top",  censor=F, xlab="Age (years)", ylab="", legend.labs=c("Small", "Medium", "Large"), legend.title="", break.time.by=.5, palette = unit.palette)
plot1b<-plot1b$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1b
plot1c<-ggsurvplot(fit.fathers, data = maturation.months, xlim=c(4,6.5), censor=F, break.time.by=.5, xlab="Age (years)", ylab="", legend.labs=c("Father absent", "Father present"), legend.title="", legend="top",palette = dad.palette)
plot1c<-plot1c$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1c
plot1d<-ggsurvplot(fit.rank.inj, data = maturation.months, xlim=c(4,7.5), censor=F, break.time.by=.5, xlab="Age (years)", ylab="", legend.labs=c("No injury", "Injury"), legend.title="", legend="top")
plot1d<-plot1d$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25)) + scale_color_brewer(palette="Dark2")
plot1d

library(patchwork)
plot1a+plot1b+plot1c + plot_annotation(tag_levels = "a")
ggsave("BE_Figure1.jpg", units="cm", width=30, height=10) 

## PLOTS FOR SUPPLEMENTARY
plotS2a<-ggsurvplot(fit.rank.3.5, data = maturation.months, xlim=c(4,6.5), legend="top",  censor=F, xlab="Age (years)", ylab="Proportion yet to mature", legend.labs=c("Low", "High"), legend.title="", break.time.by=.5, palette = rank.palette)
plotS2a<-plotS2a$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plotS2a

plotS2b<-ggsurvplot(fit.unit.conception, data = maturation.months, xlim=c(4,6.5), legend="top",  censor=F, xlab="Age (years)", ylab="", legend.labs=c("Small", "Medium", "Large"), legend.title="", break.time.by=.5, palette = unit.palette)
plotS2b<-plotS2b$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plotS2b
plotS2a+plotS2b + plot_annotation(tag_levels = "a")
ggsave("BE_FigureS2.jpg", units="cm", width=18, height=10) 

## PATH ANALYSIS
library(piecewiseSEM)
library(multcompView)
library(nlme)

## PIECEWISE

path.data$Interrupted2<-path.data$Interrupted2-1
path.data$TO.Num<-as.numeric(path.data$Takeover.Range)
gelada.SEM <- psem( #create the path model
  lm(AFR ~ Maturation.Age + UnitSize_3.5 + Conception.Rank + Interrupted2 , path.data),
  glm(Takeover.Range~UnitSize_3.5 + Dad, path.data, family="binomial"),
  lm(Maturation.Age ~ UnitSize_3.5 + Conception.Rank + Dad + Takeover.Range, path.data))

sem.sum<-summary(gelada.SEM, conserve=T, .progressBar = F)
sem.sum$coefficients[8,7]
sem.sum

## INCORPORATE UNCERTAINTY
models<-list()
P.rank<-c()
P.unit<-c()
P.dad<-c()
P.to<-c()

for(i in 1:1000) {
  sim.data<-path.data
  for(j in 1:nrow(path.data)) {
    sim.data$Uncertainty[j]<-runif(n=1, min=sim.data$DownRange[j], max=sim.data$UpRange[j])
    sim.data$AFR.Uncertainty[j]<-runif(n=1, min=sim.data$Down.AFR[j], max=sim.data$Up.AFR[j])
  }
  sim.data$Maturation.Age<-round((sim.data$Maturation.Age+sim.data$Uncertainty/30.4),1)
  sim.data$AFR<-round((sim.data$AFR+sim.data$AFR.Uncertainty/30.4),1)
  models[[i]]<- psem( #create the path model
    lm(AFR ~ Maturation.Age + UnitSize_3.5 + Conception.Rank + Interrupted2 , sim.data),
    glm(Takeover.Range~UnitSize_3.5 + Dad, sim.data, family="binomial"),
    lm(Maturation.Age ~ UnitSize_3.5 + Conception.Rank + Dad + Takeover.Range, sim.data))
  sum.sem<-summary(models[[i]], show_progress=F)
  P.rank[i]<-sum.sem$coefficients[8,7]
  P.unit[i]<-sum.sem$coefficients[7,7]
  P.dad[i]<-sum.sem$coefficients[9,7]
  P.to[i]<-sum.sem$coefficients[10,7]
}

summary(models[[11]])

sum(P.rank>0.05)/1000
sum(P.unit>0.05)/1000
sum(P.dad>0.05)/1000
sum(P.to>0.05)/1000

## PLOT
colors<-c("#0072B2","#E69F00")
AFR.plot<-ggplot(data=path.data, aes(x=Maturation.Age/12, y=AFR/12)) + geom_point(aes(color=Interrupted, shape=Interrupted), size=3) + scale_color_manual(values=colors) + 
  theme_classic() + geom_smooth(color="black", method="lm", se=T) +  
  xlab("Maturation age (years)") + ylab("Age at first birth (years)") + theme(text = element_text(size=15)) + scale_y_continuous(limits = c(5,8), breaks=seq(5,8,1)) +
  scale_x_continuous(limits = c(4,7), breaks=seq(4,7,1)) + theme(legend.position = c(0.8, 0.2)) + labs(color="Takeover?",shape="Takeover?")
AFR.plot
ggsave("BE_Figure2.jpg", units="cm", width=12, height=12) 

