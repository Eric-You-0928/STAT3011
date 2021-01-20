#select all code and run it, the result will show up!  (Rstudio)
#group two written by Li Weihao
#I will use 0.0000001 as tolerance level
tol=1e-11
bisection<-function(fx,a){
                       #a is vector
  
  repeat{
    c<-0.5*(a[1]+a[2])
    if((a[2]-c)<=tol) return(c)
    else{
      if(fx(c)==0) return(c)
      if(fx(c)*fx(a[1])<0) a[2]=c else a[1]=c
    }
  }
  
}

#in the Newton Rap method, you need to know the f`(x)  and f''(x)

NewtonRap<-function(fx,df,x)  #Give the initial x0
{
  
   repeat{
    if(abs(fx(x))<tol) return(x)
    x<-x-fx(x)/df(x)
  }
}
fx<-function(x){
  log(x+x^2)/(1+x^3)
}
dfx<-function(x){
  a<-(1+2*x)/(x+x^2)*(1+x^3)-3*x*x*log(x+x^2)
  b<-(1+x^3)^2
  a/b
}
d2fx<-function(x){
  a<-(1+x^3)/((x+x^2)^2)*((2+3*x*x+8*x^3)*(x+x^2)-(1+2*x)*(1+2*x^4+x^3+2*x))-6*x^2*(1+2*x^4+x^3+2*x)
  b<-(1+x^3)*(6*x*log(x+x^2)+(1+2*x)*3*x^2/(x+x^2))-18*x^4*log(x+x^2)
  c<-a-b
  c/((1+x^3)^3)
}

ptm <- proc.time()
xb<-bisection(dfx,c(0.1,2))
proc.time()-ptm
ptm <- proc.time()
xN<-NewtonRap(dfx,d2fx,1.0)
proc.time()-ptm
mam<-fx(xb)
x<-seq(0.01,10,0.01)
y<-fx(x)
plot(y~x,type="l",main="Plot of f(x)")
abline(h=mam,col="red")
abline(v=0.08,col="red")
abline(v=2,col="red")
cat(" root of dfx calculated by bisection method is",xb,"\n",
      "root of dfx calculated by Newton method is   ",xN,"\n",
      "maximal value of fx is                       ",mam,"\n")


