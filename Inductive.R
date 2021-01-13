library(lpSolve)

#Isotonic Separation
isotonic=function(A,a,b)
{ #A=training data set. 
  relation<-matrix(0,nrow=nrow(A),ncol=nrow(A))
  for(i in 1:nrow(A)){
    for(j in 1:nrow(A)){
      if(i==j)
      { relation[i,j]=0} #reflexive reduction
      else
      {
        if(all(A[i,-ncol(A)]>=A[j,-ncol(A)]))
          relation[i,j]=1
        else
          relation[i,j]=0
      }
    }
  }
  rTransRed=relation
  for(i in 1:nrow(relation)){
    for(j in 1:nrow(relation))
    {
      if(i!=j)
      {
        for(k in 1:nrow(relation))
        {
          if(rTransRed[i,j]==1 && rTransRed[j,k]==1)
            rTransRed[i,k]=0
        }
      }
    }
  }
  numRelation=0
  for(i in 1:nrow(rTransRed))
  {
    for(j in 1:nrow(rTransRed))
    {
      if(rTransRed[i,j]==1)
      { numRelation=numRelation +1 }
    }
  }
  coe=c()
  {
    for(i in 1:nrow(A))
      if(A[i,ncol(A)]==1)
        coe=c(coe,a*-1)
    else
      coe=c(coe,b*1)
  }
  consa1<-diag(1,nrow(A),nrow(A))
  consb1=matrix(1,nrow(A))
  k=1
  consa2=matrix(0,nrow=numRelation,ncol=nrow(A))
  for(i in 1:nrow(rTransRed))
  {
    for(j in 1:ncol(rTransRed))
    {
      if(rTransRed[i,j]==1)
      {
        consa2[k,i]=1
        consa2[k,j]=-1
        k=k+1
      }
    }
  }
  consa2=rbind2(consa2,consa1)
  
  consb2=matrix(0,nrow=nrow(consa2))
  constant=rbind(consb1,consb2)
  constraint=rbind(consa1,consa2)
  direction=c()
  for(i in 1:nrow(consa1))
  {
    direction=c(direction,"<=")
    i=i+1
  }
  for(i in 1:nrow(consa2))
  {
    direction=c(direction,">=")
    i=i+1
  }
  ans=lp(direction="min",objective.in=coe,const.mat=constraint,
         const.dir=direction,const.rhs=constant)
  
  Y_opt = round(as.vector(ans$solution))
  A_opt=cbind(A[,1:ncol(A)],Y_opt)
  A0_opt=c()
  A1_opt=c()
  for(i in 1:nrow(A_opt))
  {
    if(A_opt[i,ncol(A_opt)]==1)
      A1_opt=c(A1_opt,i)
    else
      A0_opt=c(A0_opt,i)
  }
  bpts1=c()
  bpts0=c()
  for(i in 1:length(A1_opt))
  {
    count=0
    for(j in 1:length(A1_opt))
    {
      if(rTransRed[A1_opt[i],A1_opt[j]]==0)
      { count=count+1}
    }
    if(count==length(A1_opt))
    {
      bpts1=c(bpts1,A1_opt[i])
    }
  }
  
  for(i in 1:length(A0_opt))
  {
    count=0
    for(j in 1:length(A0_opt))
    {
      if(rTransRed[A0_opt[j],A0_opt[i]]==0)
        count=count+1
    }
    if(count==length(A0_opt))
    {
      bpts0=c(bpts0,A0_opt[i])
    }
  } 
  return_list=list("dominance_relation"=relation,"TransRed"=rTransRed,"A*1"=A1_opt, "A*0"=A0_opt, "Boundary1"=bpts1,
                   "Boundary0"=bpts0,"newclass"=Y_opt)
  return(return_list)
}

#Testing: predict function takes two arguments: data to be classified in a matrix form and 
# object returned by function isotonic.
predict_IS=function(data,obj,a,b,train)
{ #data parameter is the test_data without labels.
  result=c()
  for(k in 1:nrow(data))
  {
    testdata=data[k,]
    #for class 0:
    for(i in 1:length(obj$Boundary0))
    {
      count=0
      for(j in 1:(ncol(train)-1))
      {
        if(train[obj$Boundary0[i],j]>=testdata[j])
        {
          count=count+1
        }
      }
      if(count==(ncol(train)-1))
      {
        testdata[length(testdata)+1]=0
        break
      }
    }
    if(length(testdata)!=ncol(train))
    {
      for ( i in 1:length(obj$Boundary1))
      {
        count=0
        for(j in 1:(ncol(train)-1))
        {
          if(train[obj$Boundary1[i],j]<=testdata[j])
          {
            count=count+1
          }
        }
        if(count==(ncol(train)-1))
        {
          testdata[length(testdata)+1]=1
          print(length(testdata))
          break
        }
      }
    }
    if(length(testdata)!=ncol(train))
    {dist=c()
    dist0=c()
    dist1=c()
    s=0
    for(i in 1:length(obj$Boundary0))
    {  for(j in 1:(ncol(train)-1))
    {
      if(train[obj$Boundary0[i],j]-testdata[j]>0)
      {
        s=s+(train[obj$Boundary0[i],j]-testdata[j])
      }
      else
      {
        s=s
      }
    }
      dist[i]=s
    }
    dist0=a*min(as.numeric(dist))
    dist=c()
    s=0
    for(i in 1:length(obj$Boundary1))
    {  for(j in 1:(ncol(train)-1))
    {
      if(testdata[j]-train[obj$Boundary1[i],j]>0)
      {
        s=s+(testdata[j]-train[obj$Boundary1[i],j])
      }
      else
      {
        s=s
      }
    }
      dist[i]=s
    }
    dist1=b*min(as.numeric(dist))
    
    
    if(dist1 < dist0)
    {
      testdata[length(testdata)+1]=1
    }else
    {
      testdata[length(testdata)+1]=0
    }
    }
    # print(testdata[10])
    result=c(result,testdata[length(testdata)])
  }
  result=as.numeric(result)
  return(result)
}

#Funtion to calculate P-value based on Non conformity measure (alpha)
pvalue=function(alpha)
{
  m=0
  for(i in 1:(length(alpha)-1))
  {
    if(alpha[i]>=alpha[length(alpha)])
      m=m+1
  }
  return(m/length(alpha))
}

#Split the data into training set, calibration set and testing test. Apply isotonic separation on training set.
obj<-isotonic(train.data,.2,.1)

#Calculating non conformity mearuse for Calibration data.
ctrain<-calib.data
#Measure 1
{
  alpha=c()
  for (i in 1:nrow(ctrain))
  {
    reference=train.data
    if(ctrain[i,ncol(ctrain)]==1)
    {
      dist=c()
      for(j in 1:length(obj$Boundary1))
      {
        dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
      }
      num=min(dist)
      dist=c()
      for(j in 1:length(obj$Boundary0))
      {
        dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
      }
      denum=min(dist)
      
      alpha=c(alpha,num/(num+denum))
    } else
    {
      dist=c()
      for(j in 1:length(obj$Boundary0))
      {
        
        dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
        
      }
      num=min(dist)
      dist=c()
      for(j in 1:length(obj$Boundary1))
      {
        
        dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
        
      }
      denum=min(dist)
      alpha=c(alpha,num/(num+denum))   
    }
  }
}

#Test data class 1:
pvalue1=c()
ctest=test.data1
for (i in 1:nrow(ctest))
{
  reference=newtrain.data
  {
    dist=c()
    for(j in 1:length(obj$Boundary1))
    {
      dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
    }
    num=min(dist)
    dist=c()
    for(j in 1:length(obj$Boundary0))
    {
      dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
    }
    denum=min(dist)
    alpha_test=(num/(num+denum))
  } 
  newalpha=c(alpha,alpha_test)
  pvalue1=c(pvalue1,pvalue(newalpha))
}
#Test data class 0:
pvalue0=c()
ctest=test.data0
for (i in 1:nrow(ctest))
{
  reference=newtrain.data
  {
    dist=c()
    for(j in 1:length(obj$Boundary0))
    {
      dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
    }
    num=min(dist)
    dist=c()
    for(j in 1:length(obj$Boundary1))
    {
      dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
    }
    denum=min(dist)
    alpha_test=(num/(num+denum))   
  }
  newalpha=c(alpha,alpha_test)
  pvalue0=c(pvalue0,pvalue(newalpha))
}
#To predict class based on P-Value.
pred_value=c()
for(i in 1:nrow(test.data))
{
  if(pvalue0[i]>pvalue1[i])
  {
    pred_value[i]=0
  }else
  {
    pred_value[i]=1
  }
}

#Measure 2
{
  alpha=c()
  reference=train.data
  for (i in 1:nrow(ctrain))
  {
    if(ctrain[i,ncol(ctrain)]==1)
    {
      flag=0
      for(j in 1:length(obj$Boundary1))
      {
        if(all(reference[obj$Boundary1[j],-ncol(ctrain)]<=ctrain[i,-ncol(ctrain)]))
          flag=1
        else
          flag=0
      }
      if(flag==1)
        alpha=c(alpha,0)
      if(flag==0)
      {
        for(j in 1:length(obj$Boundary0))
        {
          if(all(reference[obj$Boundary0[j],-ncol(ctrain)]>=ctrain[i,-ncol(ctrain)]))
            flag=1
          else
            flag=0 
        }
        if(flag==1)
          alpha=c(alpha,1)
        else
        {
          dist=c()
          for(j in 1:length(obj$Boundary1))
          {
            dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
          }
          num=min(dist)
          dist=c()
          for(j in 1:length(obj$Boundary0))
          {
            dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
          }
          denum=min(dist)
          alpha=c(alpha,num/(num+denum))
        }
      }
    }else 
    {
      flag=0
      for(j in 1:length(obj$Boundary0))
      {
        if(all(reference[obj$Boundary0[j],-ncol(ctrain)]>=ctrain[i,-ncol(ctrain)]))
          flag=1
        else
          flag=0
      }
      if(flag==1)
      { alpha=c(alpha,0)}
      if(flag==0)
      {
        for(j in 1:length(obj$Boundary1))
        {
          if(all(reference[obj$Boundary1[j],-ncol(ctrain)]<=ctrain[i,-ncol(ctrain)]))
            flag=1
          else
            flag=0 
        }
        if(flag==1)
          alpha=c(alpha,1)
        else 
        {
          dist=c()
          for(j in 1:length(obj$Boundary0))
          { 
            dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
          }
          num=min(dist)
          dist=c()
          for(j in 1:length(obj$Boundary1))
          {
            dist=c(dist,sum(abs(ctrain[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
          }
          denum=min(dist)
          alpha=c(alpha,num/(num+denum))
        }
      }
    }
  }
}

#Test data class 1:
pvalue1=c()
ctest=test.data1
alphatest1=c()
for (i in 1:nrow(ctest))
{
  reference=train.data
  {
    if(ctest[i,ncol(ctest)]==1)
    {
      flag=0
      for(j in 1:length(obj$Boundary1))
      {
        if(all(reference[obj$Boundary1[j],-ncol(ctrain)]<=ctest[i,-ncol(ctrain)]))
          flag=1
        else
          flag=0
      }
      if(flag==1)
        alphatest1=c(alphatest1,0)
      if(flag==0)
      {
        for(j in 1:length(obj$Boundary0))
        {
          if(all(reference[obj$Boundary0[j],-ncol(ctrain)]>=ctest[i,-ncol(ctrain)]))
            flag=1
          else
            flag=0 
        }
        if(flag==1)
          alphatest1=c(alphatest1,1)
        else
        {
          dist=c()
          for(j in 1:length(obj$Boundary1))
          {
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
          }
          num=min(dist)
          dist=c()
          for(j in 1:length(obj$Boundary0))
          {
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
          }
          denum=min(dist)
          alphatest1=c(alphatest1,num/(num+denum))
        }
      }
    }else 
    {
      flag=0
      for(j in 1:length(obj$Boundary0))
      {
        if(all(reference[obj$Boundary0[j],-ncol(ctrain)]>=ctest[i,-ncol(ctrain)]))
          flag=1
        else
          flag=0
      }
      if(flag==1)
      { alphatest1=c(alphatest1,0)}
      if(flag==0)
      {
        for(j in 1:length(obj$Boundary1))
        {
          if(all(reference[obj$Boundary1[j],-ncol(ctrain)]<=ctest[i,-ncol(ctrain)]))
            flag=1
          else
            flag=0 
        }
        if(flag==1)
          alphatest1=c(alphatest1,1)
        else 
        {
          dist=c()
          for(j in 1:length(obj$Boundary0))
          { 
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
          }
          num=min(dist)
          dist=c()
          for(j in 1:length(obj$Boundary1))
          {
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
          }
          denum=min(dist)
          alphatest1=c(alphatest1,num/(num+denum))
        }
      }
    }
  }
  newalpha=c(alpha,alphatest1)
  pvalue1=c(pvalue1,pvalue(newalpha))
}

#Test data class 0:
pvalue0=c()
ctest=test.data0
alphatest0=c()
for (i in 1:nrow(ctest))
{
  reference=train.data
  {
    if(ctest[i,ncol(ctest)]==1)
    {
      flag=0
      for(j in 1:length(obj$Boundary1))
      {
        if(all(reference[obj$Boundary1[j],-ncol(ctrain)]<=ctest[i,-ncol(ctrain)]))
          flag=1
        else
          flag=0
      }
      if(flag==1)
        alphatest0=c(alphatest0,0)
      if(flag==0)
      {
        for(j in 1:length(obj$Boundary0))
        {
          if(all(reference[obj$Boundary0[j],-ncol(ctrain)]>=ctest[i,-ncol(ctrain)]))
            flag=1
          else
            flag=0 
        }
        if(flag==1)
          alphatest0=c(alphatest0,1)
        else
        {
          dist=c()
          for(j in 1:length(obj$Boundary1))
          {
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
          }
          num=min(dist)
          dist=c()
          for(j in 1:length(obj$Boundary0))
          {
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
          }
          denum=min(dist)
          alphatest0=c(alphatest0,num/(num+denum))
        }
      }
    }else 
    {
      flag=0
      for(j in 1:length(obj$Boundary0))
      {
        if(all(reference[obj$Boundary0[j],-ncol(ctrain)]>=ctest[i,-ncol(ctrain)]))
          flag=1
        else
          flag=0
      }
      if(flag==1)
      { alphatest0=c(alphatest0,0)}
      if(flag==0)
      {
        for(j in 1:length(obj$Boundary1))
        {
          if(all(reference[obj$Boundary1[j],-ncol(ctrain)]<=ctest[i,-ncol(ctrain)]))
            flag=1
          else
            flag=0 
        }
        if(flag==1)
          alphatest0=c(alphatest0,1)
        else 
        {
          dist=c()
          for(j in 1:length(obj$Boundary0))
          { 
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary0[j],-ncol(ctrain)])))
          }
          num=min(dist)
          dist=c()
          for(j in 1:length(obj$Boundary1))
          {
            dist=c(dist,sum(abs(ctest[i,-ncol(ctrain)]-reference[obj$Boundary1[j],-ncol(ctrain)])))
          }
          denum=min(dist)
          alphatest0=c(alphatest0,num/(num+denum))
        }
      }
    }
  }
  newalpha=c(alpha,alphatest0)
  pvalue0=c(pvalue0,pvalue(newalpha))
}


#To predict class based on P-Value.
pred_value=c()
for(i in 1:nrow(test.data))
{
  if(pvalue0[i]>pvalue1[i])
  {
    pred_value[i]=0
  }else
  {
    pred_value[i]=1
  }
}
