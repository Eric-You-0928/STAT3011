#need to load library:coda (for MCMC.object),MASS (for density plot)
library(coda)
library(MASS)
#Rejection method: line 5 to line 73, select all code and run it, the result will show up!  (Rstudio)
#Gibbs sampling with M-H method:line 80 to line 180 after line 64
#select all code and run it, the result will show up!  (Rstudio)
# 
#group two written by Li Weihao
f_x<-function(x){
  sin(x*pi)*(sin(pi*x*x))^20
}
f_y<-function(y){
  sin(y*pi)*(sin(2*pi*y*y))^20
}
integrate(f_x,0,1)
integrate(f_y,0,1)
const_k<-1/(0.09988587+0.1090703)

pxy_rej<-function(x,y){
  (sin(x*pi)*(sin(pi*x*x))^20+sin(y*pi)*(sin(2*pi*y*y))^20)*const_k
}
x.points <- seq(0,1,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
# make contour plot to test whether the convergence
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- pxy_rej(x.points[i],y.points[j])
  }
}
contour(x.points,y.points,z)
#find the maximum of the pxy_rej(x,y)
x_coor<-runif(100000,0,1)
y_coor<-runif(100000,0,1)
m<-max(pxy_rej(x_coor,y_coor))
m
acc_reject<-function(n){
  x_point<-NULL
  y_point<-NULL
  cnt<-0
  repeat{
    repeat{
      u<-runif(1,0,1)
      v<-runif(1,0,1)
      if(pxy_rej(u,v)>=runif(1,0,1)*m){
        cnt<-cnt+1
        x_point<-c(x_point,u)
        y_point<-c(y_point,v)
        break
      }
    }
    if(cnt==n) break
  }
  data_matrix<-cbind(x_point,y_point)
  data_matrix
}
ptm<-proc.time()
point_collect<-acc_reject(50000)
proc.time()-ptm
point_collect<-as.data.frame(point_collect)
cor(point_collect$x_point,point_collect$y_point)


for(i in 1:50000){
  points(point_collect$x_point[i],point_collect$y_point[i])
}
cor(point_collect$x_point,point_collect$y_point)

x<-seq(0,1,0.01)
y<-x
z<-outer(x,y,pxy_rej)
den3drj<-kde2d(point_collect$x_point,point_collect$y_point)
par(mfrow=c(1,2))
persp(den3drj,theta=30,phi=15,col="lightblue",ticktype="detailed")
persp(x,y,z,theta=30,phi=15,ticktype="detailed")








#Gibbs sampling with M-H method

pxy<-function(x,y){
  if(x>1|x<0|y>1|y<0)
    return(0)
  sin(x*pi)*(sin(pi*x*x))^20+sin(y*pi)*(sin(2*pi*y*y))^20
}
x.points <- seq(0,1,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
# make contour plot to test whether the convergence
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- pxy(x.points[i],y.points[j])
  }
}
contour(x.points,y.points,z,xlab = "X", ylab = "Y")
py_given_x<-function(x,y){
  if(x>1|x<0|y>1|y<0)
    return(0)
  sin(x*pi)*(sin(pi*x*x))^20+sin(y*pi)*(sin(2*pi*y*y))^20
}
px_given_y<-function(x,y){
  if(x>1|x<0|y>1|y<0)
    return(0)
  sin(x*pi)*(sin(pi*x*x))^20+sin(y*pi)*(sin(2*pi*y*y))^20
}

sample_MH_px_given_y<-function(x,y){  #use normal as proposal distribution
  x_new<-rnorm(1,x,0.7)
  acc_rate<-min(1,px_given_y(x_new,y)/px_given_y(x,y)*dnorm(x,x_new,0.7)/dnorm(x_new,x,0.7))
  if(runif(1,0,1)<acc_rate){
    jud<<-0
    return(x_new)}
  else {
    cnt_x<<-cnt_x+1
    if(jud==1){
      cnt<<-cnt+1
      jud<<-0}
    return(x)}
}
sample_MH_py_given_x<-function(x,y){  #use normal as proposal distribution
  y_new<-rnorm(1,y,0.7)
  acc_rate<-min(1,px_given_y(x,y_new)/px_given_y(x,y)*dnorm(y,y_new,0.7)/dnorm(y_new,y,0.7))
  if(runif(1,0,1)<acc_rate)
    return(y_new)
  else {
    cnt_y<<-cnt_y+1
    jud<<-1
    return(y)}
}

MH_gibbs_cor<-function(){
  n <- 50000 #抽样个数（链的长度）  
  burn_in <- 30000 #前30000个抽样按burn-in处理 
  x<-rep(0,n)
  y<-rep(0,n)
  x[1]=runif(1,0,1)
  y[1]=runif(1,0,1)
  for(i in 1:(n-1)){
    y[i+1]=sample_MH_py_given_x(x[i],y[i])
    
    x[i+1]=sample_MH_px_given_y(x[i],y[i+1])
    
  }
  b<-burn_in+1
  x_jump<-NULL
  y_jump<-NULL
  for(i in b:n){
    x_jump<-c(x_jump,x[i])
    y_jump<-c(y_jump,y[i])
    i<-i+10
  }
 # for(i in 1:length(x_jump)){
#    points(x_jump[i],y_jump[i])
 # }
 # cor_100<<-c(cor_100,cor(x_jump,y_jump))
 list_x<-mcmc(x)
  list_y<-mcmc(y)
  return(mcmc.list(list_x,list_y))
}
jud<-0
cnt<-0
cnt_x<-0
cnt_y<-0
#cor_100<-NULL           ######################################
ptm <- proc.time()   # 把程序放在中间这里--计算运行时间    #tictoc
MMCC<-MH_gibbs_cor()
proc.time()-ptm
par(mfrow=c(1,2))
den3d<-kde2d(MMCC[[1]],MMCC[[2]])
x<-seq(0,1,0.01)
y<-x
z<-outer(x,y,pxy_rej)
persp(den3d,theta=30,phi=15,col="lightblue",ticktype="detailed")
persp(x,y,z,theta=30,phi=15,ticktype="detailed")


#finish.































#########################
## below is for convergence testing : need to run the program 100 times: troublesome

#x_d<-mcmc.list(a[[1]],a[[3]])
#y_d<-mcmc.list(a[[2]],a[[4]])
#ij=5
#while(ij < 200){
#  x_d<-lappend(x_d,a[[ij]])
#  ij<-ij+2
#}
#ij=6
#while(ij < 201){
#  y_d<-lappend(y_d,a[[ij]])
#  ij<-ij+2
#}
#x_d<-as.mcmc.list(x_d)
#y_d<-as.mcmc.list(y_d)
#gelman.diag(x_d)
#gelman.diag(y_d)
#lappend <- function(lst, obj) {
#  lst[[length(lst)+1]] <- obj
#  return(lst)
#}
#x<-seq(0,1,0.01)
#y<-x
#z<-outer(x,y,pxy_rej)
#persp(x,y,z,theta=110,phi=-5,ticktype="detailed")
#persp(x,y,z,theta=35,phi=10,col="lightblue",ticktype="detailed")
#persp(x,y,z,theta=110,phi=-5,ticktype="detailed")
##persp(den3drj,theta=110,phi=-5,col="lightblue",ticktype="detailed")








