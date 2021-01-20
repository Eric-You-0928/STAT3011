#first part: line 14 to line 250 ,select all code and run it, the result will show up!  (Rstudio)
#second part: after line 250 ,select all code and run it, the result will show up!  (Rstudio)
#group two written by Li Weihao

# genetic algorithm
#popuN 种群个体数目
#pc  交叉概率
#pm  变异概率
#maxgen    保持最大所持续的代数
#GG    终止进化的代数

# 本题染色体长度为10
length_chrom=10
GG=4000
popuN=100
pm=0.2
pc=0.95
maxgen=GG/2 # it can be modified



fitness<-function(vec){
  nom<-sum(cos(vec)^4)-2*prod(cos(vec)^2)
  i<-1:length(vec)
  den<-sqrt(sum(i*vec^2))
  abs(nom/den)
}
#染色体初始化，运用value encoding
init_chrom<-function(){
  matr<-matrix(0,popuN,length_chrom)
  for(i in 1:popuN){
    repeat{
      a<-runif(length_chrom,0,10)
      if(prod(a)>=0.75) break
    }
    matr[i,]=a
  }
  matr
}
#次轮盘赌--选择
selection<-function(matr){
  chrom_fitness<-rep(0,popuN)
  for(i in 1:popuN){
    chrom_fitness[i]<-fitness(matr[i,])
  }
  chromProb_fitness<-chrom_fitness/sum(chrom_fitness)
  index<-rep(0,popuN)
  for(i in 1:popuN){       
    pick=runif(1,0,1)
    for(j in 1:popuN){
      pick=pick-chromProb_fitness[j]
      if(pick<=0){
        index[i]=j
        break
      }
    }
  }
  matr<-matr[index,]
  matr
}
selection2<-function(matr){
  chrom_fitness<-rep(0,popuN)
  for(i in 1:popuN){
    chrom_fitness[i]<-fitness(matr[i,])
  }
  index<-rep(order(chrom_fitness)[81:100],5)
  matr<-matr[index,]
  matr
}
selection3<-function(matr){  #binary Tournament selection
  ma<-NULL
  index<-NULL
  f<-NULL
  for(i in 1:popuN){
    f<-c(f,fitness(matr[i,]))
  }
  for(i in 1:popuN){
    a1<-floor(runif(1, 1,101)) 
    a2<-floor(runif(1, 1,101))
    if(f[a1]>f[a2])
      index<-c(index,a1)
    else
      index<-c(index,a2)
  } 
  
  ma<-matr[index,]
  ma
}
#染色体交叉运算
cross<-function(matr){
  for(i in 1:popuN){
    pick_first=runif(1,0,1)
    pick_second=runif(1,0,1)
    choice_one<-ceiling(popuN*pick_first)
    choice_two<-ceiling(popuN*pick_second)
    pickS<-runif(1,0,1)
    flag<-0
    if(pickS<pc){
      #这里采用随机选择XY染色体的同一个位置，然后交叉 r*x+(1-r)*x
      while(flag==0){
      pos<-ceiling(runif(1,0,length_chrom))
      r=runif(1,0,1)
      vx<-matr[choice_one,pos]
      vy<-matr[choice_two,pos]
      matr[choice_one,pos]=r*vx+(1-r)*vy
      matr[choice_two,pos]=r*vy+(1-r)*vx
      if(prod(matr[choice_one,])>0.75&prod(matr[choice_two,])>0.75)
        flag<-1
      else{
        matr[choice_one,pos]<-vx
        matr[choice_two,pos]<-vy
      }
      
      }
    }
  }
  matr
}

cross2<-function(matr){
  for(i in 1:popuN){
    pick_first=runif(1,0,1)
    pick_second=runif(1,0,1)
    choice_one<-ceiling(popuN*pick_first)
    choice_two<-ceiling(popuN*pick_second)
    pickS<-runif(1,0,1)
    flag<-0
    
    if(pickS<pc){
      #这里采用随机选择XY染色体的同一个位置，然后交叉 r*x+(1-r)*x
      while(flag==0){
        pos<-ceiling(runif(1,0,length_chrom))
        
        vx<-matr[choice_one,pos]
        vy<-matr[choice_two,pos]
        matr[choice_one,pos]=vy
        matr[choice_two,pos]=vx
        if(prod(matr[choice_one,])>0.75&prod(matr[choice_two,])>0.75)
          
          flag<-1
        else{
          matr[choice_one,pos]<-vx
          matr[choice_two,pos]<-vy
        }
        
      }
    }
  }
  matr
}
mutation2<-function(matr){
  for(i in 1:popuN){
    pick=runif(1,0,1)
    choice<-ceiling(popuN*pick)
    r<-runif(1,0,1)
    if(r<pm){
      pos<-ceiling(runif(1,0,length_chrom))
      tem<-matr[choice,pos]  
      repeat{
        if(runif(1,0,1)>0.5)
            matr[choice,pos]<-matr[choice,pos]+0.05
        else matr[choice,pos]<-matr[choice,pos]-0.05
        
        if(0<matr[choice,pos]&matr[choice,pos]<10)
          if(prod(matr[choice,])>=0.75) break
        matr[choice,pos]=tem
      }
    }
  }
  matr
}
mutation<-function(matr){
  for(i in 1:popuN){
    pick=runif(1,0,1)
    choice<-ceiling(popuN*pick)
    r<-runif(1,0,1)
    if(r<pm){
      pos<-ceiling(runif(1,0,length_chrom))
      tem<-matr[choice,pos]  
      repeat{
        if(runif(1,0,1)>0.5)
          matr[choice,pos]<-matr[choice,pos]+0.15*runif(1,0,1)
        else matr[choice,pos]<-matr[choice,pos]-0.15*runif(1,0,1)
        
        if(0<matr[choice,pos]&matr[choice,pos]<10)
          if(prod(matr[choice,])>=0.75) break
        matr[choice,pos]=tem
      }
    }
  }
  matr
}
selection4<-function(matr){
  chrom_fitness<-rep(0,popuN)
  for(i in 1:popuN){
    chrom_fitness[i]<-fitness(matr[i,])
  }
  index<-rep(order(chrom_fitness)[21:100],5)
  matr<-matr[index,]
  matr
}
main_function<-function(){
  
  mypopu_matrix<-init_chrom()
  cat("wait for seconds...we have " ,GG," generations \n")
 # mypopu_matrix<-ss
  cur_fitness<-NULL
  for(i in 1:popuN){
    cur_fitness<-c(cur_fitness,fitness(mypopu_matrix[i,]))
  }#cur_fitness 现在是一个存着100个fitness的vector
  cur_best<-max(cur_fitness)
  keep_max<-0
  save_matrix<-mypopu_matrix
  for(j in 1:GG){
    #if(j==GG/10) pc<<-0.95
    mypopu_matrix<-selection2(mypopu_matrix)
    mypopu_matrix<-cross(mypopu_matrix)
    mypopu_matrix<-mutation(mypopu_matrix)
    cur_fitness<-NULL
    for(i in 1:popuN){
      cur_fitness<-c(cur_fitness,fitness(mypopu_matrix[i,]))
    }
    new_best<-max(cur_fitness)
    if(new_best>cur_best){ 
      cur_best<-new_best
      save_matrix<-mypopu_matrix
      keep_max<-0
    }
    if(j%%2000==0){
      cat("Wait for a minute,we already evolve ",j," generation ",GG-j," generation left\n")
    }
    keep_max<-keep_max+1
    if(keep_max>maxgen){
      cat("one individual dominate ", maxgen," generation, we stop the evoloe and maximum is: ",cur_best,"\n")
      break
    }
  }
  if(keep_max<maxgen)
      cat("after",GG," generation, the maximum value is: ",cur_best,"\n")
  #save_matrix
}
#first part: evaluate our GA: by runthe main_function() 10 times
#            with different papameters and selection methods
#we can run once for demo
main_function()



#Second part: get the maximum,based on the first part good individual

indi=c(3.1504942, 3.0826354 ,3.0267547, 3.012345 ,1.8362476, 0.3642342 ,0.3423303, 0.3635644, 0.3560176, 0.3444579)
ss<-matrix(rep(indi,100),100,10,byrow=T)
length_chrom=10
GG=5000
popuN=100
pm=0.2
pc=0.95
maxgen=100000
mutation2<-function(matr){
  for(i in 1:popuN){
    pick=runif(1,0,1)
    choice<-ceiling(popuN*pick)
    r<-runif(1,0,1)
    if(r<pm){
      pos<-ceiling(runif(1,0,length_chrom))
      tem<-matr[choice,pos]  
      repeat{
        if(runif(1,0,1)>0.5)
          matr[choice,pos]<-matr[choice,pos]+0.05*runif(1,0,1)
        else matr[choice,pos]<-matr[choice,pos]-0.05*runif(1,0,1)
        
        if(0<matr[choice,pos]&matr[choice,pos]<10)
          if(prod(matr[choice,])>=0.75) break
        matr[choice,pos]=tem
      }
    }
  }
  matr
}
main_function<-function(ss){
  # mypopu_matrix<-init_chrom()
   cat("wait for seconds...we have " ,GG," generations \n")
  mypopu_matrix<-ss
  cur_fitness<-NULL
  for(i in 1:popuN){
    cur_fitness<-c(cur_fitness,fitness(mypopu_matrix[i,]))
  }#cur_fitness 现在是一个存着100个fitness的vector
  cur_best<-max(cur_fitness)
  keep_max<-0
  save_matrix<-mypopu_matrix
  for(j in 1:GG){
    
    mypopu_matrix<-selection2(mypopu_matrix)
    mypopu_matrix<-cross(mypopu_matrix)
    mypopu_matrix<-mutation2(mypopu_matrix)
    cur_fitness<-NULL
    for(i in 1:popuN){
      cur_fitness<-c(cur_fitness,fitness(mypopu_matrix[i,]))
    }
    new_best<-max(cur_fitness)
    if(new_best>cur_best){ 
      cur_best<-new_best
      save_matrix<-mypopu_matrix
      keep_max<-0
    }
     if(j%%2000==0){
        cat("Wait for a minute,we already evolve ",j," generation ",GG-j," generation left\n")
     }
    keep_max<-keep_max+1
    if(keep_max>maxgen){
      cat("one individual dominate ", maxgen," generation, we stop the evoloe and maximum is: ",cur_best,"\n")
      break
    }
  }
  if(keep_max<maxgen)
    cat("after",GG," generation, the maximum value is: ",cur_best,"\n")
  save_matrix
}
ss<-main_function(ss)
#max: 0.7471992 




