setwd("D:\\Sekolah\\Kuliah\\Tugas Akhir\\Script R")

#######Input Data
data0=read.csv("D:\\Sekolah\\Kuliah\\Tugas Akhir\\Data\\Skripsi\\data kriminalitas5.csv", header=TRUE, sep=";")
data0
data=data0[,3:8]

data
#######Analisis Deskriptif
summary(data)
library(psych)
describe(data)

library("car")
scatterplot(KHMK ~ JPS, data = data)
scatterplot(KH ~ JPS, data = data)
scatterplot(KHMK ~ UMP, data = data)
scatterplot(KH ~ UMP, data = data)
scatterplot(KHMK ~ JP, data = data)
scatterplot(KH ~ JP, data = data)
scatterplot(KHMK ~ PKRB, data = data)
scatterplot(KH ~ PKRB, data = data)
#######Program Uji Korelasi Pearson
korelasi<-function(data)
{
  cat("\nUJI KORELASI PEARSON\n")
  cat("==============================================")
  alfa<-as.numeric(readline("\n\nInput nilai alfa : "))
  
  y1<-data[,1]
  y2<-data[,2]
  M<-length(y1)
  n<-M
  sy1y2<-sum(y1*y2)-(M*mean(y1)*mean(y2))
  sy1<-sqrt(sum(y1^2)-(M*(mean(y1))^2))
  sy2<-sqrt(sum(y2^2)-(M*(mean(y2))^2))
  korelasi<-sy1y2/(sy1*sy2)
  cat("\nkoefisien korelasi:",korelasi,"\n")
  t<-(korelasi*sqrt(M-2))/sqrt(1-(korelasi^2))
  v<-M-2
  ttabel<-qt(1-(alfa/2),v)
  cat("hipotesis:\n")
  cat("H0 : rho = 0\n")
  cat("H1 : rho ??? 0\n")
  p_value=round(2*pt(abs(t),v,lower.tail=FALSE),4)
  cat("\n========================================","\n")
  cat("nilai P-value = ",p_value,"\n")
  cat("==========================================","\n")
  cat("\nkesimpulan:\n")
  if(p_value<alfa)
  {
    cat("Tolak Ho\nada korelasi\n\n")
  }
  else
  {
    cat("Terima Ho\ntidak ada korelasi\n\n")
  }
}

#####################Program Identifikasi Kombiasi Orde Respon, Jumlah Knot,Titik Knot, dan Lambda Optimal Setiap Prediktor
mp<-function(x,eps=1e-006)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1])%*%t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)]%*%diag(1/diago)%*%t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}

trun<-function(prediktor,knot,orde)
{
  prediktor[prediktor<knot]<-knot
  b<-(prediktor-knot)^orde
  return(b)
}
datai<-data[,1:2]
data11<-data[,4]
data12<-data[,5]
data13<-data[,6]
data1<-cbind(datai,data11)
data2<-cbind(datai,data12)
data3<-cbind(datai,data13)

spline<-function(data)
{
  y1<-data[,1]
  y2<-data[,2]
  y<-c(y1,y2)
  M<-length(y1)
  n<-M
  cat("\n")
  P<-as.numeric(readline("Input Orde Maksimum : "))
  bb<-as.numeric(readline("Input Batas Bawah Lamda : "))
  ba<-as.numeric(readline("Input Batas Atas Lamda : "))
  h<-as.numeric(readline("Input Increment : "))
  k=0
  kecilGCV<-rep(0,k+2)
  kecilGCV[1]<-10^10000
  boptmaxx<-rep(0,k+1)
  pmaxx<-rep(0,k+1)
  repeat
  {
    
    k=k+1
    prediktor=data[,3]
    vl<-seq(bb,ba,h)
    nvl<-length(vl)
    GCV<-rep(0,nvl)
    MSE<-rep(0,nvl)
    n<-nrow(data)
    dataurut<-data[order(prediktor),1:3]
    x<-dataurut[,3]
    y<-c(dataurut[,1],dataurut[,2])
    z<-k+1
    r<-quantile(x,seq(0,1,by=1/z))
    lr<-length(r)
    R<-r[-lr]
    RR<-R[-1]
    cat("Jumlah Knot = ",k,"\n")
    for(i in 1:k)
    {
      cat("titik knots[",i,"] = ",RR[i],"\n")
    }
    p<-matrix(0,(P*P),2)
    p[,2]<-rep((1:P),P)
    for(i in 1:P)
    {
      c<-rep(i,P)
      p[((P*(i-1)+1):(P*i)),1]=c
    }
    pmax<-rep(0,2)
    gm<-rep(0,(n+1))
    GCVmin<-rep(0,(P*P))
    bopt<-rep(0,(P*P))
    for(m in 1:(P*P))
    {
      cat("\nORDE respon 1 :",p[m,1],"; ORDE respon 2
          :",p[m,2],"\n")
      for (r in 1:nvl)
      {
        Z1<-matrix(0,M,(p[m,1]+k))
        v11<-matrix(0,M,(p[m,1]))
        v21<-matrix(0,M,k)
        for(s in 1:(p[m,1]))
        {
          v11[,s]<-x^(s)
          v11[,(p[m,1])]<-x^(p[m,1])
        }
        for(j in 1:k)
        {
          v21[,j]<-trun(x,RR[j],p[m,1])
        }
        Z1[,1:(p[m,1]+k)]<-cbind(v11,v21)
        Z2<-matrix(0,M,(p[m,2]+k))
        v12<-matrix(0,M,(p[m,2]))
        v22<-matrix(0,M,k)
        for(s in 1:(p[m,2]))
        {
          v12[,s]<-x^(s)
          v12[,(p[m,2])]<-x^(p[m,2])
        }
        for(j in 1:k)
        {
          v22[,j]<-trun(x,RR[j],p[m,2])
        }
        Z2[,1:(p[m,2]+k)]<-cbind(v12,v22)
        ZZ<-rep(1,M)
        ZA<-cbind(ZZ,Z1)
        ZB<-cbind(ZZ,Z2)
        ZC<-matrix(0,M,((p[m,2]+k)+1))
        ZD<-matrix(0,M,((p[m,1]+k)+1))
        A<-cbind(ZA,ZC)
        B<-cbind(ZD,ZB)
        Z<-rbind(A,B)
        d1<-rep(0,(p[m,1]+1))
        d2<-rep(0,(p[m,2]+1))
        d3<-rep(1,k)
        D1<-c(d1,d3)
        D2<-c(d2,d3)
        d0<-c(D1,D2)
        D<-diag(d0)
        betatopi<-mp(t(Z)%*%Z+(M*vl[r]*D))%*%t(Z)%*%y
        ytopi<-Z%*%betatopi
        H<-Z%*%mp(t(Z)%*%Z+(n*vl[r]*D))%*%t(Z)
        MSE[r]<-(t(y-ytopi)%*%(y-ytopi))/M
        GCV[r]<-MSE[r]/(1-((1/M)*sum(diag(H))))^2
      }
      #tutup lambda
      mlamda<-cbind(vl,GCV,MSE)
      GCVmin[m]<-min(GCV)
      bopt[m]<-mlamda[mlamda[,2]==GCVmin[m],1]
      MSEE<-mlamda[mlamda[,2]==GCVmin[m],3]
      cat(" Nilai MSE = ",MSEE,"\n")
      cat(" Nilai GCV minimum = ",GCVmin[m],"\n")
      cat("Nilai Lambda Optimal saat GCV Minimum = ",bopt[m],"\n")
      betatopii<-mp(t(Z)%*%Z+(M*bopt[m]*D))%*%t(Z)%*%y
      ytopii<-Z%*%betatopii
      error<-y-ytopii
      ee<-cbind(ytopii,error)
      MSE<-(t(y-ytopii)%*%(y-ytopii))/length(y)
      JKT<-t(y-(mean(y)))%*%(y-(mean(y)))
      JKG<-t(y-ytopii)%*%(y-ytopii)
      RK<-1-(JKG/JKT)
    }
    #tutup orde
    for(a in 1:(P*P))
    {
      if(GCVmin[a]==min(GCVmin))
        
      {
        kecilGCV[k+1]<-GCVmin[a]
        boptmaxx[k]<-bopt[a]
        pmaxx[k]<-a
      }
    }
    if(kecilGCV[k+1]>kecilGCV[k])
    {
      kmax<-k-1
      boptmax<-boptmaxx[k-1]
      pmax<-pmaxx[k]
      cat("\n\nOptimal")
      cat("\nNilai GCV minimum adalah",kecilGCV[k])
      cat("\npada nilai lambda optimal = ",boptmax)
      cat("\nsaat orde respon 1 :",p[pmax,1],"\n")
      cat("dan orde respon 2 :",p[pmax,2],"\n")
      lr<-kmax+2
      z<-kmax+1
      r<-quantile(x,seq(0,1,by=1/z))
      lr<-length(r)
      R<-r[-lr]
      RR<-R[-1]
      k<-kmax
      cat("Jumlah Knot = ",kmax,"\ndengan ")
      for(i in 1:k)
      {
        cat("titik knots[",i,"] = ",RR[i],"\n")
      }
      Z1<-matrix(0,M,(p[m,1]+k))
      v11<-matrix(0,M,(p[m,1]))
      v21<-matrix(0,M,k)
      for(s in 1:(p[m,1]))
      {
        v11[,s]<-x^(s)
        v11[,(p[m,1])]<-x^(p[m,1])
      }
      for(j in 1:k)
      {
        v21[,j]<-trun(x,RR[j],p[m,1])
      }
      Z1[,1:(p[m,1]+k)]<-cbind(v11,v21)
      Z2<-matrix(0,M,(p[m,2]+k))
      v12<-matrix(0,M,(p[m,2]))
      v22<-matrix(0,M,k)
      for(s in 1:(p[m,2]))
      {
        v12[,s]<-x^(s)
        v12[,(p[m,2])]<-x^(p[m,2])
      }
      for(j in 1:k)
      {
        v22[,j]<-trun(x,RR[j],p[m,2])
      }
      Z2[,1:(p[m,2]+k)]<-cbind(v12,v22)
      ZZ<-rep(1,M)
      ZA<-cbind(ZZ,Z1)
      ZB<-cbind(ZZ,Z2)
      ZC<-matrix(0,M,((p[m,2]+k)+1))
      ZD<-matrix(0,M,((p[m,1]+k)+1))
      A<-cbind(ZA,ZC)
      B<-cbind(ZD,ZB)
      Z<-rbind(A,B)
      d1<-rep(0,(p[m,1]+1))
      d2<-rep(0,(p[m,2]+1))
      d3<-rep(1,k)
      D1<-c(d1,d3)
      D2<-c(d2,d3)
      d0<-c(D1,D2)
      D<-diag(d0)
      betatopiii<-mp(t(Z)%*%Z+(M*boptmax*D))%*%t(Z)%*%y
      ytopiii<-Z%*%betatopiii
      errorfix<-y-ytopiii
      ee<-cbind(ytopiii,errorfix)
      cat("\nPenduga Parameter : ")
      cat("\n")
      for(m in 1:((p[pmax,1]+kmax)+1))
      {
        cat("\nNilai Phi-topi [",m,"] Model 1=
            ",format(betatopiii[m]))
      }
      for(m in ((p[pmax,1]+kmax)+2):(((p[pmax,1]+kmax)+1)+
                                     ((p[pmax,2]+kmax)+1)))
      {
        cat("\nNilai Phi-topi [",m-((p[pmax,1]+kmax)+1),"]
            Model 2= ",format(betatopiii[m]))
      }
      cat("\n\n\tHasil Estimasi:\n\tytopi\t\terror\n")
      print(ee)
      cat("\n")
      MSE<-(t(y-ytopiii)%*%(y-ytopiii))/length(y)
      cat("MSE=",MSE,"\n")
      JKT<-t(y-(mean(y)))%*%(y-(mean(y)))
      JKG<-t(y-ytopiii)%*%(y-ytopiii)
      RK<-1-(JKG/JKT)
      cat("R-square=",RK,"\n")
      break
      }
    else
    {
      cat("lanjut tambah jumlah knot\n\n")
    }
    }#tutup repeat k
}

#Program Uji Glesjer dan Penentuan Matrik Pembobot
glesjer<-function(ERfix)
{
  cat("\n\nUJI GLESJER\n")
  cat("==============================================")
  alfa<-as.numeric(readline("\n\nInput nilai alfa : "))
  ###############################################
  y1<-data[,1]
  y2<-data[,2]
  y<-c(y1,y2)
  M<-length(y1)
  n<-M
  cat("\n")
  para<-as.numeric(readline("Input Banyak Prediktor Parametrik: "))
  bb<-as.numeric(readline("Input Batas Bawah Lamda : "))
  ba<-as.numeric(readline("Input Batas Atas Lamda : "))
  h<-as.numeric(readline("Input Increment : "))
  kolom=length(data[1,])
  q=ncol(data)-para-2
  k<-rep(0,q)
  orde1<-rep(0,q)
  orde2<-rep(0,q)
  z<-rep(0,q)
  lr<-rep(0,q)
  for(i in 1:q)
  {
    cat("\nInput Orde Respon 1 Prediktor ke-",i," = ")
    orde1[i]<-as.numeric(readline(" "))
    cat("Input Orde Respon 2 Prediktor ke-",i," = ")
    orde2[i]<-as.numeric(readline(" "))
    cat("Input Banyak Knot Prediktor ke-",i," = ")
    k[i]<-as.numeric(readline(" "))
    lr[i]<-k[i]+2
    z[i]<-k[i]+1
  }
  
  dataurut<-matrix(0,M,q)
  t<-matrix(0,M,q)
  RR<-matrix(0,max(k),q)
  prediktor=data[,(para+3):kolom]
  dataA=as.matrix(prediktor)
  p<-cbind(orde1,orde2)
  r<-matrix(0,max(lr),q)
  R<-matrix(0,(max(lr)-1),q)
  for(i in 1:q)
  {
    dataurut[,i]<-sort(dataA[,i])
    t[,i]<-dataurut[,i]
    if(lr[i]<max(lr))
    {
      r[1:lr[i],i]<-quantile(t[,i],seq(0,1,by=1/z[i]))
      r[(lr[i]+1):max(lr),i]<-0
    }
    else
    {
      r[,i]<-quantile(t[,i],seq(0,1,by=1/z[i]))
    }
    R[,i]<-r[-lr[i],i]
    RR[,i]<-R[-1,i]
  }
  for (i in 1:q)
  {
    cat("\nPrediktor ke-",i)
    cat("\nORDE respon 1 :",p[i,1],"; ORDE respon 2
        :",p[i,2],"\n")
    for(j in 1:k[i])
    {
      cat("titik knot [",j,"] = ",RR[j,i],"\n")
    }
  }
  y<-c(data[,1],data[,2])
  vl<-seq(bb,ba,h)
  nvl<-length(vl)
  GCV<-rep(0,nvl)
  MSE<-rep(0,nvl)
  for (r in 1:nvl)
  {
    Z1<-matrix(0,M,(sum(p[,1])+sum(k)))
    for(u in 1:q)
    {
      v11<-matrix(0,M,(p[u,1]))
      v21<-matrix(0,M,(k[u]))
      for(s in 1:p[u,1])
      {
        v11[,s]<-dataA[,u]^(s)
        v11[,(p[u,1])]<-dataA[,u]^(p[u,1])
      }
      for(j in 1:k[u])
      {
        v21[,j]<-trun(dataA[,u],RR[j,u],p[u,1])
        
      }
      if(u==1)
      {
        Z1[,1:(p[u,1]+k[u])]<-cbind(v11,v21)
      }
      else
      {
        Z1[,((sum(p[1:(u-1),1])+sum(k[1:(u-1)])+1):(sum(p[1:(u),1])+sum(k[1:(u)])))]<-cbind(v11,v21)
      }
    }
    Z2<-matrix(0,M,(sum(p[,2])+sum(k)))
    for(u in 1:q)
    {
      v12<-matrix(0,M,(p[u,2]))
      v22<-matrix(0,M,(k[u]))
      for(s in 1:p[u,2])
      {
        v12[,s]<-dataA[,u]^(s)
        v12[,(p[u,2])]<-dataA[,u]^(p[u,2])
      }
      for(j in 1:k[u])
      {
        v22[,j]<-trun(dataA[,u],RR[j,u],p[u,2])
      }
      if(u==1)
      {
        Z2[,1:(p[u,2]+k[u])]<-cbind(v12,v22)
      }
      else
      {
        Z2[,((sum(p[1:(u-1),2])+sum(k[1:(u-1)])+1):(sum(p[1:(u),2])+sum(k[1:(u)])))]<-cbind(v12,v22)
      }
    }
    ZZ<-rep(1,M)
    ZA<-cbind(ZZ,Z1)
    ZB<-cbind(ZZ,Z2)
    ZC<-matrix(0,M,((sum(p[,2])+sum(k))+1))
    ZD<-matrix(0,M,((sum(p[,1])+sum(k))+1))
    A<-cbind(ZA,ZC)
    B<-cbind(ZD,ZB)
    Z<-rbind(A,B)
    D1<-rep(0,(sum(p[,1])+sum(k)))
    D2<-rep(0,(sum(p[,2])+sum(k)))
    for(i in 1:q)
    {
      d1<-rep(0,(p[i,1]))
      for(j in 1:(p[i,1]))
      {
        d1[j]<-0
      }
      d2<-rep(0,(p[i,2]))
      for(j in 1:(p[i,2]))
      {
        d2[j]<-0
      }
      d3<-rep(0,(k[i]))
      for(j in 1:k[i])
      {
        d3[j]<-1
      }
      if(i==1)
      {
        D1[1:((p[i,1])+(k[i]))]<-c(d1,d3)
        D2[1:((p[i,2])+(k[i]))]<-c(d2,d3)
      }
      else
      {
        D1[(sum(p[1:(i-1),1])+sum(k[1:i-1])+1):(sum(p[1:i,1])+sum(k[1:i]))]<-c(d1,d3)
        D2[(sum(p[1:(i-1),2])+sum(k[1:i-1])+1):(sum(p[1:i,2])+sum(k[1:i]))]<-c(d2,d3)
      }
    }
    DD1<-c(0,D1)
    DD2<-c(0,D2)
    d0<-c(DD1,DD2)
    D<-diag(d0)
    c1<-matrix(0,M,1+para)
    for(i in 1:(para+1))
    {
      c1[,i]<-1
      c1[,1+para]<-data[,2+para]
    }
    c2<-matrix(0,M,1+para)
    for(i in 1:(para+1))
    {
      c2[,i]<-1
      c2[,1+para]<-data[,2+para]
    }
    c3<-matrix(0,M,1+para)
    c4<-matrix(0,M,1+para)
    X1<-cbind(c1,c3)
    X2<-cbind(c4,c2)
    X<-rbind(X1,X2)
    I<-diag(1,M+M)
    A<-matrix(0,M+M,M+M)
    A<-Z%*%mp(t(Z)%*%Z+(n*vl[r]*D))%*%t(Z)
    Apar<-X%*%mp(t(X)%*%t(I-A)%*%(I-A)%*%X)%*%t(X)%*%t(I-A)%*%(I-A)
    Anon<-A%*%(I-Apar)
    Asemi<-Apar+Anon
    ytopi<-Asemi%*%y
    MSE[r]<-(t(y-ytopi)%*%(y-ytopi))/n
    GCV[r]<-MSE[r]/(1-((1/n)*sum(diag(1-A))))^2
  }
  mlamda<-cbind(vl,GCV,MSE)
  GCVmin<-min(GCV)
  bopt<-mlamda[mlamda[,2]==GCVmin,1]
  MSEE<-mlamda[mlamda[,2]==GCVmin,3]
  cat("\nNilai MSE = ",MSEE,"\n")
  cat("Nilai GCV minimum = ",GCVmin,"\n")
  cat("Nilai Lambda Optimal saat GCV Minimum = ",bopt,"\n")
  I<-diag(1,M+M)
  AA<-matrix(0,M+M,M+M)
  AA<-Z%*%mp(t(Z)%*%Z+(n*bopt*D))%*%t(Z)
  AApar<-X%*%mp(t(X)%*%t(I-AA)%*%(I-AA)%*%X)%*%t(X)%*%t(I-AA)%*%(I-AA)
  AAnon<-AA%*%(I-AApar)
  AAsemi<-AApar+AAnon
  teta<-mp(X)%*%AApar%*%y
  ystar<-y-(X%*%teta)
  phi<-mp(t(Z)%*%Z+(n*bopt*D))%*%t(Z)%*%ystar
  cat("\nPenduga Parameter : ")
  for(m in 1:(para+1))
  {
    cat("\nNilai Teta-topi[",m-1,"] Model Respon 1=",format(teta[m]))
  }
  for(m in (para+2):(2*(para+1)))
  {
    cat("\nNilai Teta-topi[",m-1-(para+1),"] Model Respon 2=",format(teta[m]))
  }
  cat("\n")
  for(m in 1:((sum(p[,1])+sum(k))+1))
  {
    cat("\nNilai Phi-topi [",m-1,"] Model Respon 1=",format(phi[m]))
  }
  for(m in ((sum(p[,1])+sum(k))+2):(sum(p)+(2*sum(k))+2))
  {
    cat("\nNilai Phi-topi [",m-1-((sum(p[,1])+sum(k))+1),"]
        Model Respon 2= ",format(phi[m]))
  }
  cat("\n\n\tHasil Estimasi:\n\tytopi\t\terror\n")
  ytopii<-AAsemi%*%y
  error<-y-ytopii
  ee<-cbind(ytopii,error)
  print(ee)
  cat("\n")
  MSE<-(t(y-ytopii)%*%(y-ytopii))/length(y)
  cat("MSE=",MSE,"\n")
  JKT<-t(y-(mean(y)))%*%(y-(mean(y)))
  JKG<-t(y-ytopii)%*%(y-ytopii)
  RK<-1-(JKG/JKT)
  cat("R-square=",RK,"\n")
  ER<-matrix(0,M,2)
  ER[,1]<-error[1:M]
  ER[,2]<-error[(M+1):(2*M)]
  cat("\nNilai error untuk respon 1 dan 2 adalah \n\n")
  print(ER)
  Asemi=as.matrix(Asemi)
  error=abs(error)
  error=as.matrix(error)
  errorbar=mean(error)
  n=nrow(error)
  yhat=Asemi%*%error
  j=nrow(data)-2
  error1=error-yhat
  SSE=sum((error-yhat)^2)
  SSR=sum((yhat-errorbar)^2)
  SST=SSR+SSE
  MSE=SSE/(n-j)
  MSR=SSR/(j-1)
  cat("hipotesis:\n")
  cat("H0 : var(1)=var(2)=...=var(n)= var\nKasus
      Homoskedastisitas\n")
  cat("H1 : minimal ada satu var(i) ??? var\nKasus
      Heteroskedastisitas\n")
  Fhit=MSR/MSE
  pvalue=pf(Fhit,(j-1),(n-j),lower.tail=FALSE)
  cat("\nAnalysis of Variance","\n")
  cat("==============================================================","\n")
  cat("Sumber df SS MS Fhit
      pvalue","\n")
  cat("Regresi ",j-1," ",SSR," ",MSR,"",Fhit,"",pvalue,"\n")
  cat("Error ",n-j," ",SSE,"",MSE,"\n")
  cat("Total ",n-1," ",SST,"\n")
  cat("==============================================================","\n")
  cat("\nkesimpulan:")
  if (pvalue<=alfa)
  {
    cat("\nTolak Ho\nKasus Heteroskedastisitas","\n")
    cat("","\n")
    c<-rep(0,(M+1))
    vr1<-rep(0,M)
    vr2<-rep(0,M)
    cv<-rep(0,M)
    vr<-var(ER)
    vr1<-rep(vr[1,1],M)
    vr2<-rep(vr[2,2],M)
    cv<-rep(vr[1,2],M)
    A<-diag(vr1,M)
    B<-diag(cv,M)
    C<-B
    D<-diag(vr2,M)
    AA<-cbind(A,B)
    
    BB<-cbind(C,D)
    U<-rbind(AA,BB)
    W<-solve(U)
  }
  else
  {
    cat("\nTerima Ho\nKasus Homoskedastisitas","\n")
    cat("","\n")
    vr1<-rep(0,M)
    vr2<-rep(0,M)
    cv<-rep(0,M)
    vr<-var(ER)
    vr1<-rep(vr[1,1],M)
    vr2<-rep(vr[2,2],M)
    cv<-rep(vr[1,2],M)
    A<-diag(vr1,M)
    B<-diag(cv,M)
    C<-B
    D<-diag(vr2,M)
    AA<-cbind(A,B)
    BB<-cbind(C,D)
    U<-rbind(AA,BB)
    W<-solve(U)
  }
  print(W)
  }

##Program Estimasi Model Regresi Semiparametrik Birespon Berdasarkan Estimator Penalized Spline (dengan Pembobot)
mp<-function(x,eps=1e-006)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1])%*%t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)]%*%diag(1/diago)%*%t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}
trun<-function(prediktor,knot,orde)
{
  prediktor[prediktor<knot]<-knot
  b<-(prediktor-knot)^orde
  return(b)
}

spline2<-function(data)#dengan pembobot
{
  y1<-data[,1]
  y2<-data[,2]
  y<-c(y1,y2)
  M<-length(y1)
  n<-M
  para<-as.numeric(readline("Input Banyak Prediktor Parametrik : "))
  bb<-as.numeric(readline("Input Batas Bawah Lamda : "))
  ba<-as.numeric(readline("Input Batas Atas Lamda : "))
  h<-as.numeric(readline("Input Increment : "))
  kolom=length(data[1,])
  q=ncol(data)-para-2
  k<-rep(0,q)
  orde1<-rep(0,q)
  orde2<-rep(0,q)
  z<-rep(0,q)
  lr<-rep(0,q)
  for(i in 1:q)
  {
    cat("\nInput Orde Respon 1 Prediktor ke-",i," = ")
    orde1[i]<-as.numeric(readline(" "))
    cat("Input Orde Respon 2 Prediktor ke-",i," = ")
    orde2[i]<-as.numeric(readline(" "))
    cat("Input Banyak Knot Prediktor ke-",i," = ")
    k[i]<-as.numeric(readline(" "))
    lr[i]<-k[i]+2
    z[i]<-k[i]+1
  }
  
  dataurut<-matrix(0,M,q)
  t<-matrix(0,M,q)
  RR<-matrix(0,max(k),q)
  prediktor=data[,(para+3):kolom]
  dataA=as.matrix(prediktor)
  p<-cbind(orde1,orde2)
  r<-matrix(0,max(lr),q)
  R<-matrix(0,(max(lr)-1),q)
  for(i in 1:q)
  {
    dataurut[,i]<-sort(dataA[,i])
    t[,i]<-dataurut[,i]
    if(lr[i]<max(lr))
    {
      r[1:lr[i],i]<-quantile(t[,i],seq(0,1,by=1/z[i]))
      r[(lr[i]+1):max(lr),i]<-0
    }
    else
    {
      r[,i]<-quantile(t[,i],seq(0,1,by=1/z[i]))
    }
    R[,i]<-r[-lr[i],i]
    RR[,i]<-R[-1,i]
  }
  for (i in 1:q)
  {
    cat("\nPrediktor ke-",i)
    cat("\nORDE respon 1 :",p[i,1],"; ORDE respon 2
        :",p[i,2],"\n")
    for(j in 1:k[i])
    {
      cat("titik knot [",j,"] = ",RR[j,i],"\n")
    }
  }
  y<-c(data[,1],data[,2])
  vl<-seq(bb,ba,h)
  nvl<-length(vl)
  GCV<-rep(0,nvl)
  MSE<-rep(0,nvl)
  for (r in 1:nvl)
  {
    Z1<-matrix(0,M,(sum(p[,1])+sum(k)))
    for(u in 1:q)
    {
      v11<-matrix(0,M,(p[u,1]))
      v21<-matrix(0,M,(k[u]))
      for(s in 1:p[u,1])
      {
        v11[,s]<-dataA[,u]^(s)
        v11[,(p[u,1])]<-dataA[,u]^(p[u,1])
      }
      for(j in 1:k[u])
      {
        v21[,j]<-trun(dataA[,u],RR[j,u],p[u,1])
      }
      if(u==1)
      {
        Z1[,1:(p[u,1]+k[u])]<-cbind(v11,v21)
      }
      else
      {
        Z1[,((sum(p[1:(u-1),1])+sum(k[1:(u-1)])+1):(sum(p[1:(u),1])+sum(k[1:(u)])))]<-cbind(v11,v21)
      }
    }
    Z2<-matrix(0,M,(sum(p[,2])+sum(k)))
    for(u in 1:q)
    {
      v12<-matrix(0,M,(p[u,2]))
      v22<-matrix(0,M,(k[u]))
      for(s in 1:p[u,2])
      {
        v12[,s]<-dataA[,u]^(s)
        v12[,(p[u,2])]<-dataA[,u]^(p[u,2])
      }
      for(j in 1:k[u])
      {
        v22[,j]<-trun(dataA[,u],RR[j,u],p[u,2])
      }
      if(u==1)
      {
        Z2[,1:(p[u,2]+k[u])]<-cbind(v12,v22)
      }
      else
      {
        Z2[,((sum(p[1:(u-1),2])+sum(k[1:(u-1)])+1):(sum(p[1:(u),2])+sum(k[1:(u)])))]<-cbind(v12,v22)
      }
    }
    ZZ<-rep(1,M)
    ZA<-cbind(ZZ,Z1)
    ZB<-cbind(ZZ,Z2)
    ZC<-matrix(0,M,((sum(p[,2])+sum(k))+1))
    ZD<-matrix(0,M,((sum(p[,1])+sum(k))+1))
    A<-cbind(ZA,ZC)
    B<-cbind(ZD,ZB)
    Z<-rbind(A,B)
    D1<-rep(0,(sum(p[,1])+sum(k)))
    D2<-rep(0,(sum(p[,2])+sum(k)))
    for(i in 1:q)
    {
      d1<-rep(0,(p[i,1]))
      for(j in 1:(p[i,1]))
      {
        d1[j]<-0
      }
      d2<-rep(0,(p[i,2]))
      for(j in 1:(p[i,2]))
      {
        d2[j]<-0
      }
      d3<-rep(0,(k[i]))
      for(j in 1:k[i])
      {
        d3[j]<-1
      }
      if(i==1)
      {
        D1[1:((p[i,1])+(k[i]))]<-c(d1,d3)
        D2[1:((p[i,2])+(k[i]))]<-c(d2,d3)
      }
      else
      {
        D1[(sum(p[1:(i-1),1])+sum(k[1:i-1])+1):(sum(p[1:i,1])+sum(k[1:i]))]<-c(d1,d3)
        D2[(sum(p[1:(i-1),2])+sum(k[1:i-1])+1):(sum(p[1:i,2])+sum(k[1:i]))]<-c(d2,d3)
      }
    }
    DD1<-c(0,D1)
    DD2<-c(0,D2)
    d0<-c(DD1,DD2)
    D<-diag(d0)
    c1<-matrix(0,M,1+para)
    for(i in 1:(para+1))
    {
      c1[,i]<-1
      c1[,1+para]<-data[,2+para]
    }
    c2<-matrix(0,M,1+para)
    for(i in 1:(para+1))
    {
      c2[,i]<-1
      c2[,1+para]<-data[,2+para]
    }
    c3<-matrix(0,M,1+para)
    c4<-matrix(0,M,1+para)
    X1<-cbind(c1,c3)
    X2<-cbind(c4,c2)
    X<-rbind(X1,X2)
    I<-diag(1,M+M)
    A<-matrix(0,M+M,M+M)
    W<-diag(0.0000148,M+M,M+M)
    A<-Z%*%mp(t(Z)%*%W%*%Z+(n*vl[r]*D))%*%t(Z)%*%W
    Apar<-X%*%mp(t(X)%*%t(I-A)%*%(I-A)%*%X)%*%t(X)%*%t(I-A)%*%(I-A)
    Anon<-A%*%(I-Apar)
    Asemi<-Apar+Anon
    ytopi<-Asemi%*%y
    MSE[r]<-(t(y-ytopi)%*%(y-ytopi))/n
    GCV[r]<-MSE[r]/(1-((1/n)*sum(diag(1-A))))^2
  }
  mlamda<-cbind(vl,GCV,MSE)
  GCVmin<-min(GCV)
  bopt<-mlamda[mlamda[,2]==GCVmin,1]
  MSEE<-mlamda[mlamda[,2]==GCVmin,3]
  cat("\nNilai MSE = ",MSEE,"\n")
  cat("Nilai GCV minimum = ",GCVmin,"\n")
  cat("Nilai Lambda Optimal saat GCV Minimum = ",bopt,"\n")
  I<-diag(1,M+M)
  AA<-matrix(0,M+M,M+M)
  AA<-Z%*%mp(t(Z)%*%W%*%Z+(n*bopt*D))%*%t(Z)%*%W
  AApar<-X%*%mp(t(X)%*%t(I-AA)%*%(I-AA)%*%X)%*%t(X)%*%t(I-AA)%*%(I-AA)
  AAnon<-AA%*%(I-AApar)
  AAsemi<-AApar+AAnon
  teta<-mp(X)%*%AApar%*%y
  ystar<-y-(X%*%teta)
  phi<-mp(t(Z)%*%W%*%Z+(n*bopt*D))%*%t(Z)%*%W%*%ystar
  cat("\nPenduga Parameter : ")
  for(m in 1:(para+1))
  {
    cat("\nNilai Teta-topi[",m-1,"] Model Respon 1=",format(teta[m]))
  }
  for(m in (para+2):(2*(para+1)))
  {
    cat("\nNilai Teta-topi[",m-1-(para+1),"] Model Respon 2=
        ",format(teta[m]))
  }
  cat("\n")
  for(m in 1:((sum(p[,1])+sum(k))+1))
  {
    cat("\nNilai Phi-topi [",m-1,"] Model Respon 1=
        ",format(phi[m]))
  }
  for(m in ((sum(p[,1])+sum(k))+2):(sum(p)+(2*sum(k))+2))
  {
    cat("\nNilai Phi-topi [",m-1-((sum(p[,1])+sum(k))+1),"]
        Model Respon 2= ",format(phi[m]))
  }
  cat("\n\n\tHasil Estimasi:\n\tytopi\t\terror\n")
  ytopii<-AAsemi%*%y
  error<-y-ytopii
  ee<-cbind(ytopii,error)
  print(ee)
  #write.csv(ee,
  #file = "D:\\Sekolah\\Kuliah\\Tugas Akhir\\Data\\Skripsi\\output1.csv", row.names = FALSE)
  cat("\n")
  MSE<-(t(y-ytopii)%*%(y-ytopii))/length(y)
  cat("MSE=",MSE,"\n")
  JKT<-t(y-(mean(y)))%*%(y-(mean(y)))
  JKG<-t(y-ytopii)%*%(y-ytopii)
  JKR<-JKT-JKG
  KTR<-JKR/90
  KTG<-JKG/95
  RK<-1-(JKG/JKT)
  cat("R-square=",RK,"\n")
  Fhit<-KTR/KTG
  pvalue=pf(Fhit,90,95,lower.tail=FALSE)
  cat("Uji Serentak=",pvalue,"\n")
  Yt<-matrix(0,M,2)
  Yt[,1]<-ytopii[1:M]
  Yt[,2]<-ytopii[(M+1):(2*M)]
  ER<-matrix(0,M,2)
  ER[,1]<-error[1:M]
  ER[,2]<-error[(M+1):(2*M)]
  urut<-c(1:M)
  datagab<-cbind(Yt,data)
  #write.csv(datagab,file = "D:\\Sekolah\\Kuliah\\Tugas Akhir\\Data\\Skripsi\\Output01.csv", row.names = FALSE)
  for(r in 1:2)
  {
    win.graph()
    plot(urut,data[,r],xlim=c(min(urut),max(urut)),ylim=c(min(c(data[,r],Yt[,r])),max(c(data[,r],Yt[,r]))),xlab="Unit Observasi",ylab="Respon")
    par(new=T)
    plot(urut,Yt[,r],type="l",col="BLUE",xlim=c(min(urut),max(urut)),ylim=c(min(c(data[,r],Yt[,r])),max(c(data[,r],Yt[,r]))),xlab="Unit Observasi",ylab="Respon")
    title("Plot Estimasi",sub=paste("\n*Plot Estimasi untuk Respon ke-",r),cex.sub = 1, font.sub = 3, col.sub = "red")
  }
  }

korelasi(data)
spline(data1)
spline(data2)
spline(data3)
spline(data)
glesjer(ERfix)
spline2(data)
