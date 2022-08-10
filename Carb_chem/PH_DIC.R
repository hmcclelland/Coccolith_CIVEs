
## Set working directory to source file location. 
## Load Seacarb package.
## Run all 



# const pH : 

carbout = carb(1, 8, seq(3E-6, 500E-6, length.out = 1000), T = 25)
PH = carbout$pH
DIC = carbout$DIC
carbout1 = cbind(DIC, PH)
colnames(carbout1) = c('DIC', 'PH')
write.table(carbout1, 'carboutph.csv', sep = ',')
carbout

# const omega calcite: 
x  = seq(0.208E-3,0.250E-3,length.out = 10)

for (val in c(1:10)){
  carbout = carb(3, seq(3E-6, 500E-6, length.out = 1000), x[val], T = 25)
  PH = carbout$pH
  DIC = carbout$DIC
  carbout1 = cbind(DIC, PH)
  colnames(carbout1) = c('DIC', 'PH')
  name = paste(val,'_carboutoc.csv',sep="")
  write.table(carbout1, name, sep = ',')
  carbout
}

# const ALK: 

carbout = carb(4, seq(3E-6, 500E-6, length.out = 1000), 2E-3, T = 25)
PH = carbout$pH
DIC = carbout$DIC
carbout1 = cbind(DIC, PH)
colnames(carbout1) = c('DIC', 'PH')
write.table(carbout1, 'carboutal.csv', sep = ',')
carbout

