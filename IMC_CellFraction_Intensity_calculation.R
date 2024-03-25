# load all data
rm(list=ls())
load('IMC_LUAD.RData')

# Step 1. Cell fraction in all/immune cells, and Cell intensity calculation
myList = CellPositionList
mysiz = IMCsize
mysiz = mysiz$size_1 * mysiz$size_2/1e6

for(k in 1:length(myList)){
  cat("\r", k)
  data = myList[[k]]
  colnames(data)= c("posx", "posy", "dist.x", "dist.y", "cel.typ")
  se = which(data$dist.y<=30)
  data = data[se,]
  tmp = table(data$cel.typ)
  if(k==1){
    res = matrix(0, length(myList), length(tmp))
    row.names(res) = names(myList)
    colnames(res) = names(tmp)
    res = as.data.frame(res)
  }
  xx = tmp[colnames(res)]
  xx[is.na(xx)]=0
  res[k,] = xx
}

xx = apply(res, 1, sum, na.rm=T)
IMC_fraction = res/xx
IMC_intensity = res/mysiz
tmp = res[, c(1:2, 4:6, 8:16)]
xx = apply(tmp, 1, sum, na.rm=T)
IMC_fraction_immune = tmp/xx

colnames(IMC_intensity) = colnames(IMC_fraction) = gsub(' cell', '', gsub('MAC','Mac',gsub('Alt','M2',colnames(IMC_intensity)  )))
colnames(IMC_fraction_immune) = gsub(' cell', '', gsub('MAC','Mac',gsub('Alt','M2',colnames(IMC_fraction_immune)  )))

