# load all data
rm(list=ls())
load('IMC_LUAD.RData')

# Step 4. calculate cell count of each cell type proximal by M2 Mac with 10-100 pixels
cel.typ = c("Cancer", "Tc", "Endothelial cell", "Cl MAC", "B cell", "Th", "Alt MAC", "Treg", 
            "T other", "Neutrophils", "Cl Mo", "Mast cell", "Non-Cl Mo", "Int Mo", "NK cell", "DCs cell")
min.cell.number = 10 		## if a cell type has <10 cells, don't set it as the reference cell type in a sample
maxs = 100*10  ## the max size x 10
cel.pair.name = outer(cel.typ, cel.typ, paste, sep="__")
cel.pair.name = as.vector(t(cel.pair.name))		##  A__B --> A is reference cell type
cel.typ.num = length(cel.typ)
nnk = 1
template.mymat = matrix(0, cel.typ.num, 11)
row.names(template.mymat) = cel.typ
colnames(template.mymat) = c("All", "D10", "D20", "D30", "D40","D50","D60","D70","D80","D90","D100")

final.List = list(NULL)

for(ff in 1:length(CellPositionList)){
  cat("\r", ff)
  data = CellPositionList[[ff]]
  colnames(data) = gsub("pos\\.", "", colnames(data))
  se = which(data$dist.y<=30)
  data = data[se,]
  mylab = data$cel.typ		## the cell type label
  myx = data$posx
  myy = data$posy
  tmpx = outer(myx, myx, "-")^2
  tmpy = outer(myy, myy, "-")^2
  dist = sqrt(tmpx+tmpy)			# the distance matrix
  diag(dist) = NA
  mymat = template.mymat
  
  se = which(mylab=="Alt MAC")
  if(length(se)==0){
    xx = table(mylab)
    mymat[,1] = xx[cel.typ]
    final.List[[ff]] = mymat
    next
  }
  dist = dist[,se, drop=F]  #### important, in this, the M2s are in columns
  
  for(i in 1:cel.typ.num){
    se = which(mylab==cel.typ[i])
    if(length(se)<1){
      mymat[i,]=0
      next
    }
    dist.sub = dist[se, ,drop=F]		
    myvec = floor(apply(dist.sub, 1, min, na.rm=T))
    mymat[i,1] = nrow(dist.sub)
    mymat[i,2] = sum(myvec<=10)
    mymat[i,3] = sum(myvec<=20)
    mymat[i,4] = sum(myvec<=30)
    mymat[i,5] = sum(myvec<=40)
    mymat[i,6] = sum(myvec<=50)
    mymat[i,7] = sum(myvec<=60)
    mymat[i,8] = sum(myvec<=70)
    mymat[i,9] = sum(myvec<=80)
    mymat[i,10] = sum(myvec<=90)
    mymat[i,11] = sum(myvec<=100)
  }
  final.List[[ff]] = mymat
}
names(final.List) = names(CellPositionList)
M2ProximalCellCountList = final.List   #### M2 proximal cell counts list
rm(final.List)

