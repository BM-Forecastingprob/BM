#'Creates The Use of Marginal Distributions in Conditional Forecasting
#' @export
#' @param dt The parameter that containing the data
#' @param m number of time series
ff=function(dt=data.frame())
{
  n=nrow(dt)
  print("Enter the number of time series")
  m=scan("",quiet=TRUE)
  print("Enter the number of predicted values")
  w=scan("",quiet=TRUE)
  mx=matrix(nrow=m,ncol=1)
  mn= matrix(nrow=m,ncol=1)
  rn= matrix(nrow=m,ncol=1)
  k=2.5*(n^(0.25))
  l= matrix(nrow=m,ncol=1)
  if((round(k,0))<k)
  {k=round(k,0)+1}else
  {
    k=(round(k,0))
  }
  for(i in 1:m)
  {
    mx[i]=max(dt[,i])
    mn[i]=min(dt[,i])
    rn[i]=mx[i]-mn[i]
  }
  rn=round(rn,2)
  
  for(i in 1:m)
  {
    l[i]=(rn[i]/k)
  }
  l=round(l,2)
  xd=matrix(nrow=k,ncol=m)
  xu= matrix(nrow=k,ncol=m)
  xx=matrix(nrow=k,ncol=m)
  f= matrix(nrow=k,ncol=m)
  for(i in 1:m)
  {
    xd[1,i]=mn[i]
    xu[1,i]=mn[i]+l[i]
    xx[1,i]=(xd[1,i]+xu[1,i])/2
    for(j in 2:k)
    {
      xd[j,i]=xd[j-1,i]+l[i]
      xu[j,i]=xu[j-1,i]+l[i]
      xx[j,i]=(xd[j,i]+xu[j,i])/2
    }}
  xx=round(xx,2)
  for(i in 1:m)
  {xu[k,i]=xu[k,i]+0.5}
  for(i in 1:m)
  {
    f[,i]=0
    for(s in 1:k)
    {
      for(j in 1:n)
      {
        if (dt[j,i]>=xd[s,i] & dt[j,i]<xu[s,i])
        {f[s,i]=f[s,i]+1}
      }}}
  r=matrix(nrow=k,ncol=k^(m-1))
  rr=matrix(nrow=1,ncol=k^(m-1))
  if(m==2)
  {f2(q,dt,xd,xu,xx,k,n,r,rr,m,w)}
  if(m==3)
  {f3(q,dt,xd,xu,xx,k,n,r,rr,m,w)}
}
f2=function(q,dt,xd,xu,xx,k,n,r,rr,m,w)
{
  print("Enter independent time series values")
  q=matrix(c(scan(,,quiet=TRUE)),w,m-1)
  for(i in 1:k)
  {
    for(s in 1:k)
    {
      r[i,s]=0
      for(j in 1:n)
      {
        ss1= dt[j,1]>=xd[s,1] & dt[j,1]<xu[s,1]
        ss2= dt[j,2]>=xd[i,2] & dt[j,2]<xu[i,2]
        ss=(ss1 & ss2)
        if(ss)
        {r[i,s]=r[i,s]+1}
      }}}
  r=r/n
  r=round(r,2)
  rr=apply(r,2,sum)
  rr3=apply(r,1,sum)
  print("Common distribution table for time series"); print(r)
  print("Table of marginal distribution of independent time series");print(rr)
  print("Table of marginal distribution of dependent time series");print(rr3)
  for(i in 1:w)
  { y=0
  sss=0
  for(s in 1:k)
  {
    sss=sss+1
    if(q[i]>=xd[s,1] & q[i]<xu[s,1])
    {if(rr[sss]!=0)
    {
      y=y+(r[,sss]*xx[,2])/rr[sss]}
      else
      {for(sss in 1:k)
        y=y+(rr3[sss]*xx[sss,2])}}
    
    su=sum(y)
  }
  print("Enter the predicted values");print(su)
  }}
f3=function(q,dt,xd,xu,xx,k,n,r,rr,m,w,rr3)
{ print("Enter independent time series values")
  q=matrix(c(scan("",quiet=TRUE)),w,m-1,byrow = T)
  sss=0
  for(i in 1:k)
  {
    for(s in 1:k)
    {
      sss=sss+1
      for(h in 1:k)
      {
        r[h,sss]=0
        for(j in 1:n)
        {
          ss1= (dt[j,1]>=xd[i,1] & dt[j,1]<xu[i,1])
          ss2= (dt[j,2]>=xd[s,2] & dt[j,2]<xu[s,2])
          ss3= (dt[j,3]>=xd[h,3] & dt[j,3]<xu[h,3])
          ss=(ss1 & ss2 & ss3)
          if(ss)
          {r[h,sss]=r[h,sss]+1}
        }}}}
  r=r/n
  r=round(r,2)
  rr=apply(r,2,sum)
  rr3=apply(r,1,sum)
  print("Common distribution table for time series"); print(r)
  print("Table of marginal distribution of independent time series");print(rr)
  print("Table of marginal distribution of dependent time series");print(rr3)
  for(i in 1:w)
  {
    y=0
    sss=0
    for(s in 1:k)
    {
      for(h in 1:k)
      {
        sss=sss+1
        if(q[i,1]>=xd[s,1] & q[i,1]<xu[s,1])
        {
          if(q[i,2]>=xd[h,2] & q[i,2]<xu[h,2])
          {
            if(rr[sss]!=0)
            {
              y=y+(r[,sss]*xx[,3])/rr[sss]}
            else
            {for(sss in 1:k)
              y=y+(rr3[sss]*xx[sss,3])}}
          su=sum(y)
        }}}
    print("predicted values");print(su)
  }}

