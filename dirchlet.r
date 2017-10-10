#Sample from a Dirchlet distribution and a Dirchhlet process
require(Compositional)
require(MCMCpack)

# Sampling from a Bivriate distribution -----------------------------------
draws1= rdirichlet(800,c(1,1,1)) #uniform
draws2= rdirichlet(800,c(.1,.1,.1)) #Concentrated on tails
draws3= rdirichlet(800,c(3,3,3)) #Center, multinomial
draws4= rdirichlet(800,c(10,1,1)) #skewed
par(mfrow=c(2,2))
bivt.contour(draws1)
title("800 random draws for a D(1,1,1)")
bivt.contour(draws2)
title("800 random draws for a D(.1,.1,.1)")
bivt.contour(draws3)
title("800 random draws for a D(3,3,3)")
bivt.contour(draws4)
title("800 random draws for a D(10,1,1)")

# Building from function -------------------------------------------------------------

x = seq(0.01,0.99,length=100)
y = x
fd = function(d1,d2,astar){
  xx=cbind(d1,d2,1-d1-d2)
  return (ddirichlet(xx,alpha = astar))
}
fd(0.1,0.1,c(1,2,2))
z = outer(x,y,"fd",astar=c(1,1,1))
z1 = outer(x,y,"fd",astar=c(2,15,2))
z2 = outer(x,y,"fd",astar=c(3,3,3))
z3 = outer(x,y,"fd",astar=c(10,1,1))

par(mfrow=c(1,1))
z[is.na(z)] = 0
z1[is.na(z1)] = 0
z2[is.na(z2)] = 0
z3[is.na(z3)] = 0
op = par(bg="grey")
persp(x,y,z,theta = 50,phi = 20,expand = 0.5,zlab = "",ylab="x2",xlab="x1",col="red",main="Density of a D(1,1,1)")
persp(x,y,z1,theta = 50,phi = 20,expand = 0.5,zlab = "",ylab="x2",xlab="x1",col="red",main="Density of a D(2,15,1)")
persp(x,y,z2,theta = 50,phi = 20,expand = 0.5,zlab = "",ylab="x2",xlab="x1",col="red",main="Density of a D(3,3,3)")
persp(x,y,z3,theta = 50,phi = 20,expand = 0.5,zlab = "",ylab="x2",xlab="x1",col="red",main="Density of a D(10,1,1)")



# Sampling from a Dirchilet process ---------------------------------------

sample.dir = function(nn=10,M=5){
  x = seq(-4,4,length=11)
  y =c()
  y[1] = pnorm(x[1])
  for (i in 2:11) y[i]=pnorm(x[i])-pnorm(x[i-1])
  y=c(y,1-pnorm(x[11]))
  param = M*y
  sample.dir = rdirichlet(nn,param)
  draw= apply(t(sample.dir),2,cumsum)
  return((draw))
}
draws = sample.dir()
xx = c(seq(-4,4,length=11),5)
matplot(xx,draws,col=1:10,type = "b")
curve(pnorm(x),add = T)

draws1 = sample.dir(10,M=0.1)
matplot(xx,draws1,col=1:10,type = "b")
curve(pnorm(x),add = T)


draws2 = sample.dir(10,M=15)
matplot(xx,draws2,col=1:10,type = "b")
curve(pnorm(x),add = T)

draws3 = sample.dir(10,M=50)
matplot(xx,draws3,col=1:10,type = "b")
curve(pnorm(x),add = T)








