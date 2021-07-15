## LOAD PACKAGES NEEDED
require(plyr);require(dplyr);require(lme4); require(lmerTest); require(MuMIn); require(ggplot2);require(ggplot2);require(ggthemes);require(ggbeeswarm);require(effects);require(coxme);require(reshape2)
theme_set(theme_minimal(base_size=20))

## FUNCTION FOR ASSESSING VIFs FROM LME4 OBJECTS
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

setwd("~/Desktop/Dissertation Stuff/Jacob Pilot Fall 2018/Gelada Data")


## SURVIVAL ANALYSIS

maturation.months<-read.csv("BE_Data_1.csv")
maturation.months$End.Age[maturation.months$Start.Age==maturation.months$End.Age]<-maturation.months$End.Age[maturation.months$Start.Age==maturation.months$End.Age] + 0.001
surv_object<-Surv(time = maturation.months$Start.Age/12, time2=maturation.months$End.Age/12, event=maturation.months$Matured)

## FIT MODEL WITHOUT TIME TRANSFORMATION
fit.coxph.no.time <- coxph(surv_object ~ Takeover.Range + Rank_Conception + UnitSize_Conception + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)

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
fit.coxph1.poly <- coxph(surv_object ~ Takeover.Range + Rank_Conception + poly(UnitSize_Conception,2) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph2.poly <- coxph(surv_object ~ Takeover.Range + Rank_Conception + poly(UnitSize_3.5,2) +tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph3.poly <- coxph(surv_object ~ Takeover.Range + Rank_3.5 + poly(UnitSize_Conception,2) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph4.poly <- coxph(surv_object ~ Takeover.Range + Rank_3.5 +  poly(UnitSize_3.5,2) + tt(Dad) + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
fit.coxph.NULL <- coxph(surv_object ~ 1,  data = maturation.months)

MuMIn::model.sel(fit.coxph1.poly, fit.coxph2.poly, fit.coxph3.poly, fit.coxph4.poly, fit.coxph.NULL)
fit.coxph1.poly

## GENERATE PREDICTED CURVES
## RANK CURVES
library(broom)
library(tidyr)
library(magrittr)

fit.coxph.no.time <- coxph(surv_object ~ Takeover.Range + Rank_Conception + poly(UnitSize_Conception,2) + Dad + scale(Avg.Min) + MajorInjury + scale(Cumulative.Rainfall) + cluster(Mom_ID),  data = maturation.months)
newdf <- maturation.months %$%
  expand.grid(UnitSize_Conception=6, Rank_Conception=c(0,1), Takeover.Range=0, Dad="N", Cumulative.Rainfall=mean(maturation.months$Cumulative.Rainfall), Avg.Min=mean(maturation.months$Avg.Min), MajorInjury="N")
newdf$group <- as.character(1:nrow(newdf))

survcurv1 <- survfit(fit.coxph.no.time, newdata = newdf) %>%
  tidy() %>% 
  gather('key', 'value', -time, -n.risk, -n.event, -n.censor) %>% 
  mutate(group = substr(key, nchar(key), nchar(key)),
         key   = substr(key, 1, nchar(key) - 2)) %>% 
  left_join(newdf, 'group') %>% 
  spread(key, value)

survcurv1
ggplot(survcurv1, aes(x = time, y = estimate, ymin = conf.low, ymax = conf.high,
                     col = as.character(Rank_Conception), fill = as.character(Rank_Conception))) +
  geom_line(size = 1)+
  geom_ribbon(alpha = 0.1, col = NA) + xlim(3.5,6)

## UNIT SIZE CURVES
newdf <- maturation.months %$%
  expand.grid(UnitSize_Conception=c(2,6,10), Rank_Conception=0.5, Takeover.Range=0, Dad="N", Cumulative.Rainfall=mean(maturation.months$Cumulative.Rainfall), Avg.Min=mean(maturation.months$Avg.Min), MajorInjury="N")
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
                      col = as.character(UnitSize_Conception), fill = as.character(UnitSize_Conception))) +
  geom_line(size = 1)+
  geom_ribbon(alpha = 0.1, col = NA) + xlim(3.5,6)

## GENERATE SURVIVAL CURVES FOR PLOTS
fit.rank.conception <- survfit(surv_object ~ RankCat_Conception, data = maturation.months)
fit.rank.3.5 <- survfit(surv_object ~ RankCat_3.5, data = maturation.months)
fit.unit.conception <- survfit(surv_object ~ Unit_Cat_Conception, data = maturation.months)
fit.unit.3.5 <- survfit(surv_object ~ Unit_Cat_3.5, data = maturation.months)
fit.fathers <- survfit(surv_object ~ Dad, data = maturation.months)

rank.palette<-c("#21908CFF", "#440154FF")
unit.palette<-c("#fec44f", "#ec7014","#993404")
dad.palette<-c("#1F78B4", "#33A02C")

plot1a<-ggsurvplot(fit.rank.conception, data = maturation.months, xlim=c(4,6.5), censor=F, break.time.by=.5, xlab="Age (years)", ylab="Proportion yet to mature", legend.labs=c("Low-ranking", "High-ranking"), legend.title="", legend="top", palette = rank.palette)  
plot1a<-plot1a$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1a
plot1b<-ggsurvplot(fit.unit.conception, data = maturation.months, xlim=c(4,6.5), legend="top",  censor=F, xlab="Age (years)", ylab="", legend.labs=c("Small", "Medium", "Large"), legend.title="", break.time.by=.5, palette = unit.palette)
plot1b<-plot1b$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1b
plot1c<-ggsurvplot(fit.fathers, data = maturation.months, xlim=c(4,6.5), censor=F, break.time.by=.5, xlab="Age (years)", ylab="", legend.labs=c("Father present", "Father absent"), legend.title="", legend="top",palette = dad.palette)
plot1c<-plot1c$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plot1c

plot1a+plot1b+plot1c + plot_annotation(tag_levels = "a")

## PLOTS FOR SUPPLEMENTARY
plotS2a<-ggsurvplot(fit.rank.3.5, data = maturation.months, xlim=c(4,6.5), legend="top",  censor=F, xlab="Age (years)", ylab="", legend.labs=c("Low", "High"), legend.title="", break.time.by=.5, palette = rank.palette)
plotS2a<-plotS2a$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plotS2a

plotS2b<-ggsurvplot(fit.unit.3.5, data = maturation.months, xlim=c(4,6.5), legend="top",  censor=F, xlab="Age (years)", ylab="", legend.labs=c("Small", "Medium", "Large"), legend.title="", break.time.by=.5, palette = unit.palette)
plotS2b<-plotS2b$plot + theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))
plotS2b

## PATH ANALYSIS
path.data<-read.csv("BE_Data_2.csv")
path.data$Interrupted<-as.factor(path.data$Interrupted)
path.data$Interrupted2<-as.numeric(path.data$Interrupted)

library(piecewiseSEM)
library(multcompView)
help(package="piecewiseSEM")
gelada.SEM <- psem( #create the path model
  lm(AFR ~ Maturation.Age + Conception.Group + Conception.Rank + Interrupted2, path.data),
  lm(Maturation.Age ~ Conception.Group + Conception.Rank + Dad + Takeover.Range, path.data))

summary(gelada.SEM, conserve=T)
AIC(gelada.SEM)
plot(gelada.SEM)

colors<-c("#0072B2","#E69F00")
AFR.plot<-ggplot(data=path.data, aes(x=Maturation.Age/12, y=AFR/12)) + geom_point(aes(color=Interrupted, shape=Interrupted), size=3) + scale_color_manual(values=colors) + 
  theme_classic() + geom_smooth(color="black", method="lm", se=T) + 
  xlab("Maturation age (years)") + ylab("Age at first birth (years)") + theme(text = element_text(size=15)) + scale_y_continuous(limits = c(5,8), breaks=seq(5,8,1)) +
  scale_x_continuous(limits = c(4,6.2), breaks=seq(3,7,1)) + theme(legend.position = c(0.85, 0.2)) + labs(color="Takeover?",shape="Takeover?")
AFR.plot
