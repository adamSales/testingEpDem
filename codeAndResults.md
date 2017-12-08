---
title: Replication Code and Results for "Testing Epistemic Democracy's Claims for Majority Rule"
author: Adam C Sales and Zev Berger
format: html_document
---




```r
library(lme4)
library(stargazer)
library(splines)
#print(load('datNonImputed.RData'))
```



```r
if(!is.element('anes',ls())){
anes <- read.table('data/anes_timeseries_cdf_rawdata.txt', sep='|',header=TRUE)
}

## keep only years of interest
anes <- anes[anes$VCF0004%in%seq(1956,2012,4),]

anes$vote <- anes$VCF0706
anes$vote[anes$vote%in%c(0,7)] <- NA

dat <- data.frame(year=anes$VCF0004)

dat$weight <- anes$VCF0009x

dat$state <- anes$VCF0901b

dat$repub <- dat$year %in% c(1956,1960,1972,1976,1984,1988,1992,2004,2008)

dat$inc <- ifelse(dat$repub,2,1)

dat$vote <- anes$vote==dat$inc


dat$vote <- rep(NA,nrow(dat))
dat$vote[dat$repub & anes$vote==2] <- 1
dat$vote[dat$repub & anes$vote==1] <- 0
dat$vote[!dat$repub & anes$vote==2] <- 0
dat$vote[!dat$repub & anes$vote==1] <- 1



dat$incumbentParty <- ifelse(dat$repub,-1,1)
dat$finances <- ifelse(anes$VCF0880==1,1,
                       ifelse(anes$VCF0880==2,0,
                              ifelse(anes$VCF0880==3,-1,NA)))
dat$finances1 <- ifelse(anes$VCF0880b==1,1,
                       ifelse(anes$VCF0880b==3,0,
                              ifelse(anes$VCF0880b==5,-1,NA)))
dat$finances[dat$year<1966] <- dat$finances1[dat$year<1966]

dat$partyID <- ifelse(anes$VCF0301==1,3,
                      ifelse(anes$VCF0301%in%c(2,3),2,
                             ifelse(anes$VCF0301%in%c(4,9),0, #Apolitical=ind??
                                    ifelse(anes$VCF0301%in%c(5,6),-2,
                                           ifelse(anes$VCF0301==7,-3,NA)))))
dat$region <- anes$VCF0112

dat$retrospective <- anes$VCF0871
dat$retrospective <- ifelse(dat$retrospective%in%c(8,9),
                            ifelse(anes$VCF0870==1, 2,
                                   ifelse(anes$VCF0870==5,4,NA)),dat$retrospective)
dat$retrospective[dat$retrospective==0] <- NA
dat$retrospective <- -(dat$retrospective-3)/2

dat$prospective <- -(anes$VCF0872-3)/2
dat$prospective[dat$prospective%in%c(-2.5,-3)] <- NA

dat$race <- as.numeric(anes$VCF0106!=1)
dat$race[anes$VCF0106==9] <- NA
dat$race6 <- factor(anes$VCF0105b)
dat$race6[dat$race6==0] <- 9
#dat$race6[dat$race6%in%c('3','4')] <- '7'



nbiRdi <- read.csv('data/NBI.RDI.csv') #this is cut and pasted from Table 1 in Nadeau Lewis-Beck.
## I'm assuming the data they put in the table is the data they used in the analysis.

dat <- merge(dat,nbiRdi,by='year',all.x=TRUE)

dat$finances1 <- NULL
dat$repub <- NULL


### add in some more variables to help with the imputation
dat$age <- anes$VCF0101
dat$age <- dat$age/sd(dat$age,na.rm=TRUE)
dat$gender <- factor(anes$VCF0104,exclude=NULL)
dat$hisp <- factor(anes$VCF0107  )
dat$hisp[dat$hisp%in%c('1','2','3','4')] <- '1'
dat$educ <- factor(anes$VCF0110,ordered=TRUE,exclude=NULL)
dat$urbanism <- factor(anes$VCF0111,exclude=NULL)
dat$marital <- anes$VCF0147
dat$marital[dat$marital==4] <- 3
dat$marital[dat$marital==7] <- 2
dat$marital[dat$marital>7] <- NA
dat$marital <- factor(dat$marital)
dat$class <- factor(anes$VCF0148,exclude=NULL)
dat$class[dat$class=='7'] <- '9'



yearExt=c(2000,2004,2008,2012)
rdiDat <- read.csv('data/RDI.csv',na.st='.',stringsAsFactors=FALSE)
# downloaded from FRED http://research.stlouisfed.org/fred2/series/A229RX0/downloaddata 5/27/14


nbiDat <- read.csv('data/NBIquarterly.csv',skip=1)
# downloaded from http://www.sca.isr.umich.edu/data-archive/mine.php 'Table 25: Current Business Conditions Compared with a Year Ago' 5/27/14
nbiDat$index <- nbiDat[,'Better.Now']-nbiDat[,'Worse.Now']


nbiYear <- function(year){

    row <- function(quarter)
        which(nbiDat$Year==year & nbiDat$Quarter==quarter)

    nbiDat$index[row(4)]
}

nbiRdiExtend <- cbind(
    year=yearExt,
    NBI=sapply(yearExt, nbiYear),

    RDI=sapply(yearExt,function(year)
        rdiDat[rdiDat$DATE==paste(year,'-11-01',sep=''),'VALUE'])
)

for(y in c(2000,2004,2008,2012)){
    row <- nbiRdiExtend[nbiRdiExtend[,'year']==y,]
    dat$NBI[dat$year==y] <- row['NBI']
    dat$RDI[dat$year==y] <- row['RDI']
}
```

## Models reported in paper


```r
rf <- glmer(vote~RDI+finances+incumbentParty+incumbentParty:race+incumbentParty:partyID+(1|year)+(1|state),
            data=dat,family=binomial(logit),weights=weight)
ret <- update(rf,.~.+retrospective)
pro <- update(rf,.~.+prospective)
both <- update(rf,.~.+retrospective+prospective)
```


```r
stargazer(rf,ret,pro,both,type='html',star.cutoffs=c(0.05,0.01,0.001),
          notes=' <sup> * </sup>p<0.05; <sup> ** </sup>p<0.01; <sup> *** </sup>p<0.001',notes.append=FALSE)
```


<table style="text-align:center"><tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="4"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="4" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="4">vote</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td><td>(3)</td><td>(4)</td></tr>
<tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">RDI</td><td>0.207<sup>***</sup></td><td>0.025</td><td>0.200<sup>***</sup></td><td>0.028</td></tr>
<tr><td style="text-align:left"></td><td>(0.060)</td><td>(0.053)</td><td>(0.057)</td><td>(0.057)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">finances</td><td>0.305<sup>***</sup></td><td>0.180<sup>***</sup></td><td>0.316<sup>***</sup></td><td>0.185<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.030)</td><td>(0.041)</td><td>(0.042)</td><td>(0.044)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">incumbentParty</td><td>-0.490<sup>***</sup></td><td>-0.466<sup>***</sup></td><td>-0.325<sup>***</sup></td><td>-0.454<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.102)</td><td>(0.087)</td><td>(0.094)</td><td>(0.093)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">retrospective</td><td></td><td>1.126<sup>***</sup></td><td></td><td>1.066<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td>(0.075)</td><td></td><td>(0.080)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">prospective</td><td></td><td></td><td>0.435<sup>***</sup></td><td>0.293<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td>(0.050)</td><td>(0.052)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">incumbentParty:race</td><td>1.188<sup>***</sup></td><td>1.082<sup>***</sup></td><td>1.056<sup>***</sup></td><td>1.045<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.075)</td><td>(0.091)</td><td>(0.094)</td><td>(0.096)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">incumbentParty:partyID</td><td>0.854<sup>***</sup></td><td>0.882<sup>***</sup></td><td>0.888<sup>***</sup></td><td>0.865<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.012)</td><td>(0.016)</td><td>(0.017)</td><td>(0.017)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>-0.355<sup>*</sup></td><td>0.169</td><td>-0.436<sup>**</sup></td><td>0.142</td></tr>
<tr><td style="text-align:left"></td><td>(0.179)</td><td>(0.141)</td><td>(0.148)</td><td>(0.150)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td><td></td><td></td></tr>
<tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>19,358</td><td>13,190</td><td>12,040</td><td>11,977</td></tr>
<tr><td style="text-align:left">Log Likelihood</td><td>-6,223.000</td><td>-3,189.000</td><td>-2,944.000</td><td>-2,835.000</td></tr>
<tr><td style="text-align:left">Akaike Inf. Crit.</td><td>12,462.000</td><td>6,397.000</td><td>5,907.000</td><td>5,690.000</td></tr>
<tr><td style="text-align:left">Bayesian Inf. Crit.</td><td>12,525.000</td><td>6,464.000</td><td>5,973.000</td><td>5,763.000</td></tr>
<tr><td colspan="5" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="4" style="text-align:right"><sup> * </sup>p<0.05; <sup> ** </sup>p<0.01; <sup> *** </sup>p<0.001</td></tr>
</table>

## Are there any "anti-RDI" sub-populations?


```r
dat <- within(dat,{
    financesFac <- factor(finances,exclude=NULL)
    partyIDFac <- factor(partyID,exclude=NULL)})
fullInteraction <- glmer(vote~RDI*(financesFac+partyIDFac+
                                   race6+age+class+educ+gender+marital+urbanism)+
                             incumbentParty*(race6+partyID)+
                             (RDI|state)+(1|year),data=dat,family=binomial(logit))


slopes <- fixef(fullInteraction)
slopes <- slopes[grep('RDI',names(slopes))]
slopes[-1] <- slopes[-1]+slopes['RDI']

### for whom is the coefficient on RDI negative?
terms <- names(fixef(fullInteraction))
intTerms <- terms[grep('RDI:',terms)]
intWith <- gsub('RDI:','',intTerms)

RDIslope <- slopes['RDI']+model.matrix(fullInteraction)[,intWith]%*%slopes[intTerms]+
    ranef(fullInteraction)$state[as.character(fullInteraction@frame$state),'RDI']
hist(RDIslope,main='Individual RDI Coefficients by Covariates',xlab='RDI Coefficient')
```

![plot of chunk hetero1](figure/hetero1-1.png)

## Are there any sub-populations that can't judge RDI?


```r
fullInteractionRet <- lmer(retrospective~RDI*(financesFac+partyIDFac+
                                   race6+age+class+educ+gender+marital+urbanism)+
                             incumbentParty*(race6+partyID)+
                             (RDI|state)+(1|year),data=dat)

slopes <- fixef(fullInteractionRet)
slopes <- slopes[grep('RDI',names(slopes))]
slopes[-1] <- slopes[-1]+slopes['RDI']

### for whom is the coefficient on RDI negative?
terms <- names(fixef(fullInteractionRet))
intTerms <- terms[grep('RDI:',terms)]
intWith <- gsub('RDI:','',intTerms)

RDIslope <- slopes['RDI']+model.matrix(fullInteractionRet)[,intWith]%*%slopes[intTerms]+
        ranef(fullInteractionRet)$state[as.character(fullInteractionRet@frame$state),'RDI']
hist(RDIslope,main='Individual RDI Coefficients by Covariates',xlab='RDI Coefficient')
```

![plot of chunk hetero2](figure/hetero2-1.png)

### Full interaction model results


```r
stargazer(fullInteraction,fullInteractionRet,type='html',star.cutoffs=c(0.05,0.01,0.001),
          notes=' <sup> * </sup>p<0.05; <sup> ** </sup>p<0.01; <sup> *** </sup>p<0.001',notes.append=FALSE)
```


<table style="text-align:center"><tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td>vote</td><td>retrospective</td></tr>
<tr><td style="text-align:left"></td><td><em>generalized linear</em></td><td><em>linear</em></td></tr>
<tr><td style="text-align:left"></td><td><em>mixed-effects</em></td><td><em>mixed-effects</em></td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">RDI</td><td>0.128</td><td>0.118</td></tr>
<tr><td style="text-align:left"></td><td>(0.408)</td><td>(0.203)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">financesFac0</td><td>0.293<sup>**</sup></td><td>0.157<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.113)</td><td>(0.015)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">financesFac1</td><td>0.650<sup>***</sup></td><td>0.210<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.107)</td><td>(0.013)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">financesFacNA</td><td>-0.179</td><td>0.149<sup>*</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.610)</td><td>(0.060)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">partyIDFac-2</td><td>-0.620<sup>**</sup></td><td>0.019</td></tr>
<tr><td style="text-align:left"></td><td>(0.194)</td><td>(0.021)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">partyIDFac0</td><td>-0.623<sup>**</sup></td><td>-0.009</td></tr>
<tr><td style="text-align:left"></td><td>(0.207)</td><td>(0.024)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">partyIDFac2</td><td>-0.324</td><td>0.003</td></tr>
<tr><td style="text-align:left"></td><td>(0.186)</td><td>(0.020)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">partyIDFac3</td><td>-0.237</td><td>0.039</td></tr>
<tr><td style="text-align:left"></td><td>(0.203)</td><td>(0.023)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race62</td><td>0.244</td><td>0.057<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.208)</td><td>(0.018)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race63</td><td>-0.147</td><td>0.023</td></tr>
<tr><td style="text-align:left"></td><td>(0.199)</td><td>(0.022)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race64</td><td>-0.458</td><td>0.060</td></tr>
<tr><td style="text-align:left"></td><td>(0.368)</td><td>(0.040)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race69</td><td>-0.065</td><td>0.027</td></tr>
<tr><td style="text-align:left"></td><td>(0.906)</td><td>(0.104)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">age</td><td>0.110<sup>*</sup></td><td>0.005</td></tr>
<tr><td style="text-align:left"></td><td>(0.050)</td><td>(0.006)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">class1</td><td>0.107</td><td>-0.013</td></tr>
<tr><td style="text-align:left"></td><td>(0.633)</td><td>(0.069)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">class2</td><td>-0.042</td><td>0.118</td></tr>
<tr><td style="text-align:left"></td><td>(0.823)</td><td>(0.090)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">class3</td><td>0.173</td><td>-0.014</td></tr>
<tr><td style="text-align:left"></td><td>(0.652)</td><td>(0.071)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">class4</td><td>-0.051</td><td>-0.031</td></tr>
<tr><td style="text-align:left"></td><td>(0.633)</td><td>(0.069)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">class5</td><td>-1.361</td><td>0.106</td></tr>
<tr><td style="text-align:left"></td><td>(0.828)</td><td>(0.092)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">class6</td><td>-0.143</td><td>0.007</td></tr>
<tr><td style="text-align:left"></td><td>(0.651)</td><td>(0.071)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">class9</td><td>-0.074</td><td>0.039</td></tr>
<tr><td style="text-align:left"></td><td>(0.669)</td><td>(0.072)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">classNA</td><td>0.649</td><td>0.207</td></tr>
<tr><td style="text-align:left"></td><td>(0.735)</td><td>(0.260)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">educ.L</td><td>-1.313<sup>**</sup></td><td>-0.082</td></tr>
<tr><td style="text-align:left"></td><td>(0.416)</td><td>(0.058)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">educ.Q</td><td>1.151<sup>**</sup></td><td>0.034</td></tr>
<tr><td style="text-align:left"></td><td>(0.352)</td><td>(0.049)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">educ.C</td><td>-0.755<sup>***</sup></td><td>-0.005</td></tr>
<tr><td style="text-align:left"></td><td>(0.229)</td><td>(0.032)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">educ4</td><td>0.339<sup>**</sup></td><td>0.008</td></tr>
<tr><td style="text-align:left"></td><td>(0.117)</td><td>(0.017)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">gender2</td><td>0.037</td><td>-0.016</td></tr>
<tr><td style="text-align:left"></td><td>(0.088)</td><td>(0.012)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">marital2</td><td>-0.005</td><td>0.024</td></tr>
<tr><td style="text-align:left"></td><td>(0.130)</td><td>(0.016)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">marital3</td><td>-0.074</td><td>0.007</td></tr>
<tr><td style="text-align:left"></td><td>(0.141)</td><td>(0.017)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">marital5</td><td>-0.142</td><td>0.007</td></tr>
<tr><td style="text-align:left"></td><td>(0.151)</td><td>(0.019)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">urbanism1</td><td>1.061</td><td>0.003</td></tr>
<tr><td style="text-align:left"></td><td>(1.087)</td><td>(0.769)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">urbanism2</td><td>1.050</td><td>0.008</td></tr>
<tr><td style="text-align:left"></td><td>(1.084)</td><td>(0.769)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">urbanism3</td><td>1.107</td><td>0.002</td></tr>
<tr><td style="text-align:left"></td><td>(1.085)</td><td>(0.769)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">urbanismNA</td><td>0.716</td><td>0.077</td></tr>
<tr><td style="text-align:left"></td><td>(0.651)</td><td>(0.444)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">incumbentParty</td><td>-0.620<sup>***</sup></td><td>0.087</td></tr>
<tr><td style="text-align:left"></td><td>(0.121)</td><td>(0.106)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:financesFac0</td><td>0.064</td><td>0.011<sup>*</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.040)</td><td>(0.006)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:financesFac1</td><td>0.036</td><td>0.039<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.039)</td><td>(0.005)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:financesFacNA</td><td>0.169</td><td>0.009</td></tr>
<tr><td style="text-align:left"></td><td>(0.193)</td><td>(0.019)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:partyIDFac-2</td><td>0.229<sup>***</sup></td><td>-0.009</td></tr>
<tr><td style="text-align:left"></td><td>(0.067)</td><td>(0.008)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:partyIDFac0</td><td>0.279<sup>***</sup></td><td>-0.015</td></tr>
<tr><td style="text-align:left"></td><td>(0.073)</td><td>(0.009)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:partyIDFac2</td><td>0.192<sup>**</sup></td><td>0.008</td></tr>
<tr><td style="text-align:left"></td><td>(0.065)</td><td>(0.008)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:partyIDFac3</td><td>0.194<sup>**</sup></td><td>0.005</td></tr>
<tr><td style="text-align:left"></td><td>(0.070)</td><td>(0.009)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:race62</td><td>0.045</td><td>-0.030<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.074)</td><td>(0.007)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:race63</td><td>0.079</td><td>-0.018</td></tr>
<tr><td style="text-align:left"></td><td>(0.079)</td><td>(0.009)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:race64</td><td>0.144</td><td>-0.018</td></tr>
<tr><td style="text-align:left"></td><td>(0.141)</td><td>(0.016)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:race69</td><td>0.051</td><td>-0.009</td></tr>
<tr><td style="text-align:left"></td><td>(0.295)</td><td>(0.042)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:age</td><td>-0.036<sup>*</sup></td><td>0.0003</td></tr>
<tr><td style="text-align:left"></td><td>(0.018)</td><td>(0.003)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:class1</td><td>0.062</td><td>-0.013</td></tr>
<tr><td style="text-align:left"></td><td>(0.276)</td><td>(0.034)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:class2</td><td>0.144</td><td>-0.050</td></tr>
<tr><td style="text-align:left"></td><td>(0.343)</td><td>(0.041)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:class3</td><td>0.052</td><td>-0.011</td></tr>
<tr><td style="text-align:left"></td><td>(0.281)</td><td>(0.035)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:class4</td><td>0.105</td><td>-0.004</td></tr>
<tr><td style="text-align:left"></td><td>(0.276)</td><td>(0.034)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:class5</td><td>0.674</td><td>-0.046</td></tr>
<tr><td style="text-align:left"></td><td>(0.345)</td><td>(0.042)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:class6</td><td>0.176</td><td>-0.017</td></tr>
<tr><td style="text-align:left"></td><td>(0.281)</td><td>(0.035)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:class9</td><td>0.172</td><td>-0.031</td></tr>
<tr><td style="text-align:left"></td><td>(0.288)</td><td>(0.036)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:classNA</td><td>0.091</td><td></td></tr>
<tr><td style="text-align:left"></td><td>(0.305)</td><td></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:educ.L</td><td>0.482<sup>**</sup></td><td>0.041</td></tr>
<tr><td style="text-align:left"></td><td>(0.159)</td><td>(0.023)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:educ.Q</td><td>-0.438<sup>**</sup></td><td>-0.027</td></tr>
<tr><td style="text-align:left"></td><td>(0.134)</td><td>(0.019)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:educ.C</td><td>0.269<sup>**</sup></td><td>0.009</td></tr>
<tr><td style="text-align:left"></td><td>(0.085)</td><td>(0.012)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:educ4</td><td>-0.119<sup>**</sup></td><td>-0.010</td></tr>
<tr><td style="text-align:left"></td><td>(0.042)</td><td>(0.006)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:gender2</td><td>-0.007</td><td>-0.001</td></tr>
<tr><td style="text-align:left"></td><td>(0.031)</td><td>(0.004)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:marital2</td><td>-0.025</td><td>-0.003</td></tr>
<tr><td style="text-align:left"></td><td>(0.047)</td><td>(0.006)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:marital3</td><td>0.011</td><td>0.002</td></tr>
<tr><td style="text-align:left"></td><td>(0.051)</td><td>(0.007)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:marital5</td><td>0.041</td><td>-0.003</td></tr>
<tr><td style="text-align:left"></td><td>(0.052)</td><td>(0.008)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:urbanism1</td><td>-0.295</td><td>0.025</td></tr>
<tr><td style="text-align:left"></td><td>(0.312)</td><td>(0.221)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:urbanism2</td><td>-0.343</td><td>0.019</td></tr>
<tr><td style="text-align:left"></td><td>(0.312)</td><td>(0.221)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">RDI:urbanism3</td><td>-0.360</td><td>0.027</td></tr>
<tr><td style="text-align:left"></td><td>(0.312)</td><td>(0.221)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race62:incumbentParty</td><td>2.312<sup>***</sup></td><td>0.033<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.121)</td><td>(0.010)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race63:incumbentParty</td><td>0.564<sup>***</sup></td><td>0.003</td></tr>
<tr><td style="text-align:left"></td><td>(0.093)</td><td>(0.010)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race64:incumbentParty</td><td>0.489<sup>**</sup></td><td>-0.017</td></tr>
<tr><td style="text-align:left"></td><td>(0.153)</td><td>(0.018)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">race69:incumbentParty</td><td>0.778</td><td>-0.031</td></tr>
<tr><td style="text-align:left"></td><td>(0.463)</td><td>(0.049)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">incumbentParty:partyID</td><td>0.869<sup>***</sup></td><td>0.064<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td>(0.012)</td><td>(0.002)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>-1.358</td><td>-0.686</td></tr>
<tr><td style="text-align:left"></td><td>(1.246)</td><td>(0.738)</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>20,449</td><td>20,809</td></tr>
<tr><td style="text-align:left">Log Likelihood</td><td>-7,023.000</td><td>-13,206.000</td></tr>
<tr><td style="text-align:left">Akaike Inf. Crit.</td><td>14,196.000</td><td>26,562.000</td></tr>
<tr><td style="text-align:left">Bayesian Inf. Crit.</td><td>14,791.000</td><td>27,158.000</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup> * </sup>p<0.05; <sup> ** </sup>p<0.01; <sup> *** </sup>p<0.001</td></tr>
</table>

## Which presidential popular votes would flip...
### ... were RDI at its mean value?

```r
## state-by-state too hard; look at popular vote
pop <- read.csv('data/popularVote.csv',header = TRUE)[1:15,]

### what if each year RDI had been at its mean?
ddd <- rf@frame
meanRDI <- mean(aggregate(ddd$RDI,list(year=ddd$year),mean)$x)
factual <- predict(rf)
counterfactual <- factual-fixef(rf)['RDI']*(ddd$RDI-meanRDI)
predPropFac <- aggregate(plogis(factual),by=list(year=ddd$year),mean)
predPropCouFac <- aggregate(plogis(counterfactual),by=list(year=ddd$year),mean)
pop$diff <- 100*(predPropCouFac$x-predPropFac$x)

print(with(pop,year[(-1)^(incWin)*diff> Margin/2]))
```

```
## [1] 1960 1980
```

### ... were RDI to shift by 1 SD? (in the "wrong" direction)


```r
rdisd <- sd(vapply(unique(dat$year),function(y) mean(dat$RDI[dat$year==y]),1))

pop$logodds <- log((pop$incVote/100)/(1-pop$incVote/100))
pop$counter <- ilogit(pop$logodds+(-1)^(pop$incWin)*fixef(rf)['RDI'])


print(pop$year[abs(pop$incVote/100-pop$counter)>pop$Margin/200])
```

```
##  [1] 1960 1968 1976 1980 1988 1992 1996 2000 2004 2008 2012
```


### ... were retro at its mean value?

```r
pop2 <- subset(pop,year>1979)
ddd <- ret@frame
meanRET <- mean(aggregate(ddd$retrospective,list(year=ddd$year),mean)$x)
factual <- predict(ret)
counterfactual <- factual-fixef(ret)['retrospective']*(ddd$retrospective-meanRET)
predPropFac <- aggregate(plogis(factual),by=list(year=ddd$year),mean)
predPropCouFac <- aggregate(plogis(counterfactual),by=list(year=ddd$year),mean)
pop2$diff <- 100*(predPropCouFac$x-predPropFac$x)

print(with(pop2,year[(-1)^(incWin)*diff> Margin/2]))
```

```
## [1] 1980 1992 2008
```

### ... were retro to shift by 1 SD? (in the "wrong" direction)


```r
retsd <- sd(vapply(unique(ddd$year),function(y) mean(ddd$retrospective[ddd$year==y]),1))
ilogit <- function(x) exp(x)/(1+exp(x))

pop2$logodds <- log((pop2$incVote/100)/(1-pop2$incVote/100))
pop2$counter <- ilogit(pop2$logodds+(-1)^(pop2$incWin)*fixef(ret)['retrospective'])


print(pop2$year[abs(pop2$incVote/100-pop2$counter)>pop2$Margin/200])
```

```
## [1] 1980 1984 1988 1992 1996 2000 2004 2008 2012
```




## upper bound on proportion anti-growth


```r
n <- with(dat,sum(!is.na(RDI) & !is.na(retrospective) & RDI<1.7))
succsesses <- with(dat,sum(!is.na(RDI) & !is.na(retrospective) & RDI<1.7& retrospective==1))
print(binom.test(61,5089,alternative='less',conf.level=0.95))
```

```
## 
## 	Exact binomial test
## 
## data:  61 and 5089
## number of successes = 61, number of trials = 5100, p-value
## <0.0000000000000002
## alternative hypothesis: true probability of success is less than 0.5
## 95 percent confidence interval:
##  0.00000 0.01481
## sample estimates:
## probability of success 
##                0.01199
```

```r
print(binom.test(61,5089,alternative='less',conf.level=0.99))
```

```
## 
## 	Exact binomial test
## 
## data:  61 and 5089
## number of successes = 61, number of trials = 5100, p-value
## <0.0000000000000002
## alternative hypothesis: true probability of success is less than 0.5
## 99 percent confidence interval:
##  0.00000 0.01604
## sample estimates:
## probability of success 
##                0.01199
```

```r
print(binom.test(61,5089,alternative='less',conf.level=0.999))
```

```
## 
## 	Exact binomial test
## 
## data:  61 and 5089
## number of successes = 61, number of trials = 5100, p-value
## <0.0000000000000002
## alternative hypothesis: true probability of success is less than 0.5
## 99.9 percent confidence interval:
##  0.00000 0.01748
## sample estimates:
## probability of success 
##                0.01199
```


