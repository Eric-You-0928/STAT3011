#select all code and run it, the result will show up!  (Rstudio)
#group two written by Li Weihao

quad_tri<-function(fx,a,b){   # trapezoidal rules
  m=1e7
  h=(b-a)/m
  coord<-fx(seq(a,b,by=h))
  a<-sum(coord[2:m])+coord[1]+coord[m+1]
  a*h
}

quad_sim<-function(fx,a,b){    #Simpson's rule
  m=1000000#0
  h=(b-a)/m
  coord<-seq(a,b,by=h)
  a<-sum(fx(coord[1:(m-1)])+4*fx(coord[2:m])+fx(coord[3:(m+1)]))
  a*h/6
}

fx<-function(x){abs(cos(x))/x*exp(-(log(x)-3)^2)}
x<-seq(0.01,200,0.01)
y<-fx(x)
plot(y~x,type="l")
abline(v=200,col="red")
# monte carlo integration: log-normal: mean 3  sd:sqrt(0.5)
ptm <- proc.time()   # 把程序放在中间这里--计算运行时间    #tictoc


s<-rlnorm(1e7,3,sqrt(0.5))
a<-abs(cos(s))*sqrt(pi)
proc.time()-ptm
mean(a)
ptm <- proc.time()
quad_tri(fx,0.001,200)     #quadrature method1: before 200 we find the integral
proc.time()-ptm
# done
#quadrature method2 每个奇点分开
ptm <- proc.time()
summm<-quad_tri(fx,0.001,pi/2)

for(i in 1:65){
  summm<-quad_tri(fx,pi*i-pi/2,pi*i+pi/2)+summm
} # this is time consuming: about two minutes
summm
proc.time()-ptm
f_tran<-function(t){
  
  fx(t/(1-t))/(1-t)^2
}
t<-seq(0.601,0.999,0.001)
y<-f_tran(t)
plot(y~t,type="l")
ptm<-proc.time()   #quadrature method3:transform the integral
quad_tri(f_tran,0.8,0.999)+quad_tri(f_tran,0.001,0.8)
proc.time()-ptm


d_x<-function(x){abs(cos(x))}
x<-seq(0,150,0.01)
y<-d_x(x)
plot(y~x,type="l" )

