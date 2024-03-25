rm(list=ls())
library(survival)
library(survminer)
library(ggplot2)
load('IMC_LUAD.RData')
distance_th = 30

info = IMC_clinical
cel.typ = c("Cancer", "Tc", "Endothelial cell", "Cl MAC", "B cell", "Th", "Alt MAC", "Treg", "T other", 
            "Neutrophils", "Cl Mo", "Mast cell", "Non-Cl Mo", "Int Mo", "NK cell", "DCs cell")
myList = M2ProximalCellCountList

mycol=(distance_th+10)/10		## count the cells without any M2 macrophages with (mycol-1)*100 distances  # 5 = D40
for(k in 1:length(myList))
{
  cat("\r",k)
  tmp = myList[[k]]
  if(k==1)
  {
    data = matrix(0, length(myList), nrow(tmp))
    colnames(data) = row.names(tmp)
    row.names(data) = names(myList)
    data = as.data.frame(data)
    data2 = data
  }
  all.cell.num = sum(tmp[,1], na.rm=T)
  imm.cell.num = all.cell.num - tmp["Cancer",1] - tmp["Endothelial cell",1]	
  data[k,] = (tmp[,1]-tmp[,mycol])/imm.cell.num 
  data2[k,] = tmp[,mycol]/imm.cell.num 
}

Opt = 1   # 1: all stage, 2: early, 3: late
if(Opt==1) {}
if(Opt==2) 
{
  se = which(info$Stage.Late1==0)
  info = info[se,]
}
if(Opt==3) 
{
  se = which(info$Stage.Late1==1)
  info = info[se,]
}

comxx = intersect(row.names(info), row.names(data))
data = data[comxx,]
data2 = data2[comxx,]
data[is.na(data)] = 0
data2[is.na(data2)] = 0
info = info[comxx,]
colnames(data2) = colnames(data) = gsub(' cell','',gsub('Alt ','M2 ',gsub('MAC','Mac',colnames(data))))
data = subset(data, select = -c(Cancer,Endothelial))
data2 = subset(data2, select = -c(Cancer,Endothelial))
dim(data)
dim(data2)
survreg.pval1 = coxph.pval1 = coxph.pval2 = rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 = rep(0, ncol(data))
for(k in 1:ncol(data)){
  cat("\r", k)
  mytf = as.numeric(data[,k])
  mytf2 = as.numeric(data2[,k])
  xx = cbind(mytf,mytf2, info)
  xx = xx[xx[, "time"]>0,]
  mycox = survreg(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)$table
  survreg.pval1[k] = mycox["mytf", "p"]
  mycox = coxph(Surv(time, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  hr1[k] = mycox$conf.int[1]
  lb1[k] = mycox$conf.int[3]
  ub1[k] = mycox$conf.int[4]
  mycox = coxph(Surv(time, event)~mytf2, xx) 
  mycox = summary(mycox)
  coxph.pval2[k] = mycox$coefficients[5]
  hr2[k] = mycox$conf.int[1]
}
survreg.qval1 = p.adjust(survreg.pval1, "BH")
coxph.qval1 = p.adjust(coxph.pval1, "BH")
coxph.qval2 = p.adjust(coxph.pval2, "BH")
name = colnames(data)
res1 = data.frame(name, P=coxph.pval1,  HR= hr1, QV=coxph.qval1) # without M2 approach
res2 = data.frame(name, P=coxph.pval2,  HR= hr2, QV=coxph.qval2) # with M2 approach

M2distal_survival = res1
M2proximal_survival = res2