# load all data
rm(list=ls())
load('IMC_LUAD.RData')

# Step 2. X -> Y proximity calculation and permutation
cel.typ = c("Cancer", "Tc", "Endothelial cell", "Cl MAC", "B cell", "Th", "Alt MAC", "Treg", "T other", "Neutrophils", "Cl Mo", "Mast cell", "Non-Cl Mo", "Int Mo", "NK cell", "DCs cell")
min.cell.number = 10 		## if a cell type has <10 cells, don't set it as the reference cell type in a sample
maxs = 100*10  ## the max size x 10
cel.pair.name = outer(cel.typ, cel.typ, paste, sep="__")
cel.pair.name = as.vector(t(cel.pair.name))		##  A__B --> A is reference cell type
cel.typ.num = length(cel.typ)
nnk = 1
template.mymat = matrix(0, cel.typ.num*cel.typ.num, 101)
row.names(template.mymat) = cel.pair.name
colnames(template.mymat) = c("CCI", paste("R", 1:100, sep=""))	

myList = CellPositionList

final.List = list(NULL)

for(ff in 1:length(myList))
{
  cat("\r", ff)
  data = myList[[ff]]
  colnames(data) = gsub("pos\\.", "", colnames(data))
  se = which(data$dist.y<=30)
  data = data[se,]
  mylab = data$cel.typ
  xx = table(mylab)
  mynum = xx[cel.typ]
  mynum[is.na(mynum)] = 0
  
  myx = data$posx
  myy = data$posy
  tmpx = outer(myx, myx, "-")^2
  tmpy = outer(myy, myy, "-")^2
  dist = sqrt(tmpx+tmpy)
  diag(dist) = NA
  mymat = template.mymat
  
  mycol = 1
  for(k in 1:cel.typ.num)
  {
    mysta = cel.typ.num*(k-1)+1
    myend = cel.typ.num*k
    if(mynum[k]<min.cell.number)
    {
      mymat[mysta:myend,] = NA
      next
    }
    se = which(mylab==cel.typ[k])
    dist.sub1 = dist[se,]
    for(i in 1:cel.typ.num)
    {
      if(mynum[i]<min.cell.number)
      {
        mymat[mysta+i-1,] = NA
        next
      }
      se = which(mylab==cel.typ[i])
      dist.sub2 = dist.sub1[, se]
      
      myvec = floor(apply(dist.sub2, 1, min, na.rm=T)*10)
      myvec = myvec[myvec<=maxs]
      tmp = table(myvec)
      xx = rep(0, maxs+1)
      xx[as.integer(names(tmp))+1]=tmp
      
      for(j in 2:(maxs+1))
      {
        xx[j] = xx[j-1]+xx[j]
      }
      mymat[mysta+i-1, mycol] = sum(xx/max(xx))/(maxs+1)
    }
  }
  ## 100 permutations ()
  for(p in 1:100)
  {
    cat("\r\r\r", ff, "-->", p)
    pm.mylab = sample(mylab)
    mycol = p+1
    for(k in 1:cel.typ.num)
    {
      if(mynum[k]<min.cell.number)
      {
        next
      }
      mysta = cel.typ.num*(k-1)+1
      myend = cel.typ.num*k
      se = which(pm.mylab==cel.typ[k])
      dist.sub1 = dist[se,]
      for(i in 1:cel.typ.num)
      {
        if(mynum[i]<min.cell.number)
        {
          next
        }
        se = which(pm.mylab==cel.typ[i])
        dist.sub2 = dist.sub1[, se]
        myvec = floor(apply(dist.sub2, 1, min, na.rm=T)*10)
        myvec = myvec[myvec<=maxs]
        tmp = table(myvec)
        xx = rep(0, maxs+1)
        xx[as.integer(names(tmp))+1]=tmp
        
        for(j in 2:(maxs+1))
        {
          xx[j] = xx[j-1]+xx[j]
        }
        mymat[mysta+i-1, mycol] = sum(xx/max(xx))/(maxs+1)
      }
    }	
  }
  final.List[[ff]] = mymat
}

names(final.List) = names(myList)
CCI.List = final.List  ## Proximity and permutation 100 times of all patients
rm(final.List)


# Step 3. get Proximity table of all IMC samples
df = NULL
for(i in 1:length(CCI.List)){
  tmp = CCI.List[[i]]
  df = cbind(df, tmp[,1])
}
colnames(df) = names(CCI.List)
df = as.data.frame(t(df))
se = gsub(' cell', '', gsub('MAC','Mac',gsub('Alt','M2',colnames(df)  )))
se = strsplit(se, '__')
numerator = sapply(se,"[[",2)
denominator = sapply(se,"[[",1)
se = paste0(numerator,'→',denominator)
colnames(df) = se
df = df[,!grepl('NK|DCs|Int Mo',colnames(df))]
IMC_proximity = df  #### Proximity scores
rm(df,tmp,se,numerator,denominator,i)
