source('preprocess.R')

canc.present.1 = add.vals(canc, 1, 2, 7)

canc$only1.present = rep('Present', length(canc$Gene))
for(i in seq(1,length(canc$Gene), 4))
{
  if(canc[i+3, 4] == canc[i+3,5])
  {
    canc[i,7] = 'Not detected'
    canc[i+1,7] = 'Not detected'
    canc[i+2,7] = 'Not detected'
    canc[i+3,7] = 'Not detected'
  }
}

