###Bioscience3 figures by Species
###opening stuff
setwd("G:\Quagga\Qmodel")
library(polynom)

#################
#Psuedo Code
#1) Based on a given number of life stages
#2) Create matrix model of population growth, including harvest parameters
#3) Find where the sum of the population structure yeilds a stable age distribution where its R = Harvest

#################Populatio Model structure and paramters are mostly comming from CASAGRANDI et al 2007 Zebra mussel model
###These initial parameters yeild a stable population like in CASAGRANDI Fig2e
class<-c(1,2,3,4)#a 4 stage life history, skipping velegers
N.init<-N.init<-c(0,50,50,50)#dnorm(class, mean = 2, sd = .5, log=F)*30000  #initial class distribution, follows a log-normal distribution to approximate Cope 
sig<-c(0.00001,.88,.41,.35,.04)#Natural survivorship of each age class, approx means of what CASAGRANDI considered
fert<-c(0,.24e6,.456e6,.795e6)#Fertility of each age class, CASAGRANDI fig 2
B.ini<-.001 # is filtration rate of velegers by adualts, but it's just a scaler, so doesn't matter except to get the population sizes to look right
V.init<-c(0,0,0,0) #Vulnerability to harvest
wgt<-c(0,1e-04,1e-03,1e-02) #the average weight (kg) of each class

###Market Parameters and functions 
a.Q<-200 	#the max willingness to pay for one kg
b.Q<-4500	#the max quantity in kg
c.Q<-5 	# cost of a unit of effort, on average, rather than some upstart cost
q.Q<-1e-3   #the fraction of the population harvested by one unit of effort 
M.init<-1e6 #Max value of ecosystem services
d<-0.3 #discount rate
harv<-function (a=a.Q,b=b.Q,c=c.Q,q=q.Q,N) {n<-sum(N)
	if (n<=0) return (0) 
	H<-b*(1-(c/(a*q*sum(n))))
	if (H<0) return (0) else return (H)} 
price<-function (a,b,xh) {if (xh<=0) return(a)  ##Price function: simple linear model, how much is a harvest worth to SELL?
	price<-a*(1-xh/b) 
	if (price<0) return (0) else return (price)}

###Population growth model equations
n1<-function(N=N.init,s=sig,f=fert,B=B.ini) s[1] *exp(-B*sum(N))* sum(f*N/2)       #Survival & fertility for velegers, that is the number of velegers that settle
n2<-function(s=sig,N=N.init,V=V.init,a1,b1,q1) as.numeric(s[2]*N[1]*(1-V[1]*harv(N=N,a=a1,b=b1,q=q1)/sum(N)))
n3<-function(s=sig,N=N.init,V=V.init,a1,b1,q1) as.numeric(s[3]*N[2]*(1-V[2]*harv(N=N,a=a1,b=b1,q=q1)/sum(N)))
n4<-function(s=sig,N=N.init,V=V.init,a1,b1,q1) as.numeric(s[4]*N[3]*(1-V[3]*harv(N=N,a=a1,b=b1,q=q1)/sum(N)) +s[5]*N[4]*(1-V[4]*harv(N=N,a=a1,b=b1,q=q1)/sum(N)))
#####
#####Impact Curves
yoko<-function (n,u,b,Y,M){ #these are named close to how Yokomizo names varibles
	B<-1/(1+exp(u/b))
	C<-(1+exp(-(1-u)/b)) /(1-B*(1+exp(-(1-u)/b)))
	M-M*C*(1/(1+exp(-(n/Y-u)/b))-B) }
	##M is the max value of the ecosystem service, Y is carrying Capacity, n is population, u & b are shape parameters
YY<-list(i=c(u=0,b=.1),ii=c(u=0.5,b=.1),iii=c(u=1,b=1),iv=c(u=1,b=0.1))##parameters for 4 Yoko functional forms
K.init<-50000#Carrying capacity, only for Yoko impact curves density/square meter, a large guess based on densities in Cope

##A model to simulate the populations Over time

quagga<-function(pars,N=N.init,T=25,c1=c.Q,V=V.init,y=YY[[1]],M=M.init,K=K.init,Q=q.Q,C=c.Q){    #T=time, par<-c(a=a.Q,b=b.Q), the demand parameters
	pop.t<-data.frame(n.1=N.init[1],n.2=N.init[2],n.3=N.init[3],n.4=N.init[4],P=0,H=0,Ct=0,ES=0)  #,p.2=0,p.3=0,p.4=0) 
	for (t in 2:T) { 
		Ntm1<-pop.t[t-1,1:4]
		Ht<-harv(a=pars[1],b=pars[2],c=c1,q=Q,sum(Ntm1)) #Harvest after last reproduction
		Pt<-if (Ht<=0) 0 else  price(a=pars[1],b=pars[2],Ht)
		Ct<-C/(Q*sum(V*Ntm1)) #cost of a unit of effort to harvest of the vulnerable population, if harvesting occurs.
		Nt<- c(  n1(N=pop.t[t-1,1:4]),  n2(N=pop.t[t-1,1:4],V=V,a1=pars[1],b1=pars[2],q1=Q), 
			 n3(N=pop.t[t-1,1:4],V=V,a1=pars[1],b1=pars[2],q1=Q),   n4(N=pop.t[t-1,1:4],V=V,,a1=pars[1],b1=pars[2],q1=Q)   )  
		ES<-yoko(n=sum(Nt),u=y[1],b=y[2],Y=K,M=M)#ecosystem service value at Nt
		pop.t[t,]<-c(Nt,Pt,Ht,Ct,ES)
	}
	return(pop.t)
}#end function
matplot(quagga(pars=c(a.Q,b.Q),V=c(0,.6,.6,.6))[1:4],type="b")
#####
#We want to find the DEMAND (a and b) that maximizes P*H-C*H+amenity over the time

	sim<-quagga(pars=c(10000,5000),V=c(0,.6,.6,.6))[-1,]
	value.t<-(sim$P * sim$H) - (sim$Ct*sim$H) #The undiscounted value
	sum(exp(-d*1:dim(sim)[1])*value.t)##the net present value under given conditions.
	
	










#######################
###some functions
######################
price<-function (a,b,xh) {if (xh<=0) return(a)  ##Price function: how much is a harvest worth to SELL?
	price<-a*(1-xh/b) 
	if (price<0) return (0) else return (price)}
harv<-function (a,b,c,q,xb) {if (xb<=0) return (0)  ##The Cost of harvesting a given amount
	H<-b*(1-(c/(a*q*xb)))
	if (H<0) return (0) else return (H)}  
Gr<-function (xb,r=r1,K=K1)if (xb<=0) return(0) else r*xb*(1-(xb/K))
Poly<-function (a,b,c,q,r,K) {   #this is a function to find the REAL roots of the Cubic in Eq.6, that is the Population size where harvest=growth
	roots<-polyroot(c( -((b*c)/(a*q)), b, -r, r/K   ))  #the polynomial which solves harvest=growth
	min(c(K,  max(Re(roots[abs(Im(roots))<1e-1]),na.rm=T)))#makes sure the root is less than K
	} 
yoko<-function (n,u,b,Y,M){ #these are named close to how Yokomizo names varibles
	B<-1/(1+exp(u/b))
	C<-(1+exp(-(1-u)/b)) /(1-B*(1+exp(-(1-u)/b)))
	M-M*C*(1/(1+exp(-(n/Y-u)/b))-B) }
	##M is the max value of the ecosystem service, Y is carrying Capacity, n is population, u & b are shape parameters
YY<-list(i=c(u=0,b=.1),ii=c(u=0.5,b=.1),iii=c(u=1,b=1),iv=c(u=1,b=0.1))##parameters for 4 Yoko functional forms

#a set of functions to find Xb, the population size where the net social value of population reduction+harvest = cost of harvest
social.est<-function (t,Xb,a,r,K,b,c,q,d,y,m) {# The solution makes the function=0 
	Yoko<-yoko(n=Xb,u=YY[[y]][1],b=YY[[y]][2],Y=K,M=m) 	# the environmental cost of a population of size Xb
	GofX<-r*Xb*(1-(Xb/K))						#growth of a population of sizr Xb
	CofX<-c/(q*Xb)							#the cost of a catch of size c from population of size Xb
	Pofh<-demand.prod(a=a,b=b,xh= Gr(xb=Xb,r=r,K=K)) #Price received for given harvest, assume equilibrium where harv=growth
	(exp(1)^(-d*t)) * ((Pofh*GofX - CofX*GofX )+Yoko)} #finds the long-term revenue
social.int<-function (Xb,a,r,K,b,c,q,d,y,m) 
	integrate(social.est, Xb=Xb,a=a,r=r,K=K,b=b,c=c,q=q,d=d,y=y,m=m,lower=0,upper=Inf,subdivisions=1000)$value
	#takes the integral of the long-term revenue, that is, the net present value (I think)
social.opt<-function (a,r,K,b,c,q,d,y,m) 
	optimize(social.int,interval=c(0,K),maximum=T,a=a,r=r,K=K,b=b,c=c,q=q,d=d,y=y,m=m)
	#finds the population size that maximizes the net present value
#

soc.sub<-function (m=a1*K1,a=a1,b=b1,c=c1,q=q1,r=r1,K=K1,y=2,d=d1){
	Xb<-social.opt(m=m,a=a,b=b,c=c,q=q,r=r,K=K,y=y,d=d)$maximum
	Xh<-Gr(xb=Xb,r=r,K=K)
	ch<-(c/(q*Xb))*Xh
	ph<-a*(1-Xh/b)*Xh
	esv<-yoko(n=Xb,u=YY[[y]][1],b=YY[[y]][2],Y=K,M=m)
	-m+esv+(ph-ch)
	}
d1<-.03 #discount value	
###############
#Parameters
###############
################
####crayfish, from Brett Peter's disx
################
a.cf<-4*2#1  #the max willingness to pay for one unit, lake tahoe price*2
b.cf<-450	#the max quantity, 10X more than estimated from phone call
c.cf<-320/150  #16man-hours at $20/hr
r.cf<-log(1.2)
q.cf<-2.938144e-05
k.cf<-58200  
m.cf<-k.cf*a.cf	#Relative to the maximum possible value of harvesting all individuals at one go

##cf growth function
grow.cf<-sapply(seq(1,k.cf,len=200),Gr,r=r.cf,K=k.cf)
harv.cf<-sapply(seq(1,k.cf,len=200), harv,a=a.cf,b=b.cf,c=c.cf,q=q.cf)

dev.new(width=4, height=5)
par(mfrow=c(2,1),oma=c(3,1,1,1),mar=c(1,4,0,1))
plot(y=grow.cf,x=seq(0,1,length.out = length(grow.cf)),ylab="Growth, harvest",type="l",lty=1,lwd=2,cex.lab=1,col="green",xaxt="n",xlab="")
	points(y=harv.cf,x=seq(0,1,length.out = length(grow.cf)),type="l",lwd=2,col="magenta")
	abline(v=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf)/k.cf,lty=3,col="black",lwd=1)
plot(y=yoko(n=seq(1,k.cf,len=200),u=0,b=.1,Y=k.cf,M=m.cf)/m.cf,x=seq(0,1,length.out = length(grow.cf)),lty=1, type="l",lwd=2,main="",cex.lab=1,xlab="Population Biomass",col="saddlebrown",
	ylab="Ecosystem service value",cex=.8)
	abline(v=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf)/k.cf,lwd=1,lty=3,col="black")
	arrows(x1=c(-1,-1,-1),x0=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf)/k.cf,
		y0=yoko(n=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf),u=0,b=.1,Y=k.cf,M=m.cf)/m.cf,lty=3,length=0,col="black",lwd=1)

################
####Asian Carp, silver carp, entirely from garvey et al 
################
a.ac<-.35*.453*2  #10-50cents a pound is what people are saying, so, *2 for max	
b.ac<-15000000*.453  ##15mil pound per year into meal http://www.thetelegraph.com/news/local/article_076b55f4-7bd4-11e2-9950-0019bb30f31a.html  
	(30000000/3)*.453	#to china over three years, so there's large potential demand, even if it's higher than biomass, 
	#opportunity to explore price increases
c.ac <-(20*(50+35)+10*250+10*1000+10*3000)/10   #10 fishermen and a crew member was conducted during fall 2011, 
	#resulting in the removal of about 200,000lbs
	#20*(50+35)+10*250 licences ets   +
	#10*1000 fuel
	#10*3000 gear
	#the unit of effort then is one boat
r.ac<- log(2.03) #unfished lambda, garvey et al ch12
q.ac<-(200000*.453)/(3100000*.453*.9)/10## the amount caught in the incentive program/total pop at that time
k.ac<-3100000*.453*.9 # million pounds to kg, 90% silver, Peoria Lock to the Mississippi River, garvey et al 
m.ac<-7000000000 # for the great lakes fisheries  max ESvalue 

dev.new(width=4, height=5)
par(mfrow=c(2,1),oma=c(3,1,1,1),mar=c(1,4,0,1))
grow.ac<-Gr(xb=seq(1,k.ac,length.out=200),r=r.ac,K=k.ac)
harv.ac<-sapply(seq(1,k.ac,len=200),FUN=harv,a=a.ac,b=b.ac,c=c.ac,q=q.ac)
plot(y=grow.ac,x=seq(0,1,length.out = length(grow.ac)),ylab="Growth, harvest",type="l",lty=1,lwd=2,cex.lab=1,col="green",
		xaxt="n",xlab="",ylim=c(0,max(c(harv.ac,grow.ac))))
	points(y=harv.ac,x=seq(0,1,length.out = length(harv.ac)),type="l",lwd=2,col="magenta")
	abline(v=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac)/k.ac,lty=3,col="black",lwd=1)
plot(y=yoko(n=seq(1,k.ac,length.out=200),u=1,b=0.1,Y=k.ac,M=m.ac)/m.ac,x=seq(0,1,length.out = length(grow.ac)),lty=1, type="l",lwd=2,main="",cex.lab=1,xlab="Population Biomass",col="saddlebrown",
	ylab="Ecosystem service value",cex=.8)
	abline(v=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac)/k.ac,lwd=1,lty=3,col="black")
	arrows(x1=c(-1,-1,-1),x0=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac)/k.ac,
		y0=yoko(n=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac),u=1,b=.1,Y=k.ac,M=m.ac)/m.ac,lty=3,length=0,col="black",lwd=1)

################
####Lionfish, morris et al 2011, barbour et al and others
################
a.lf<-12*.453*2  #$12/lb http://www.coralmagazine-us.com/content/jamaicans-invite-lionfish-dinner
b.lf<-20000	 #max quant, 10x more than a big derby catch.	
c.lf<-15	#Ashley an hour/dive $10 for two tanks $5for 1 hour boat fuel
		#marions$200/hour for a contract diver...assume 1 hour per dive with prep etc...cost of effort: dive team per trip
r.lf<-log(1.12) #the lambda from Morris et al 2011
q.lf<-19/(350*1240)#%K per effort:  19LF/14min with 350 lionfish per hectare, assume the surveyed the whole ha
wgt<-c(0.03,0.13,0.25,0.36,	0.45,0.52,0.56,0.59,	0.62,0.63,0.64,0.65,	0.65,0.65,0.65,0.66,	0.66,0.66,0.66,0.66)
pop<-c(100,61,37,22,14,8,5,3,2,1,1,0,0,0,0,0,0,0,0,0)##from Barbour et al 2011
k.lf<-(1240*100)*200  *sum(wgt)/sum(pop)     
		#Carrying capacty  coral reef area of 1,240 km2 corral reefs of Jamaca, but more like 200f/ha
		#reef area in ha *fish/ha *kg/fish (from barbour, average size across the size distrubtion)
m.lf<-11187937  #max ES value USD Moonsammy, S. et al. GCFI:64 (2012)

dev.new(width=4, height=5)
par(mfrow=c(2,1),oma=c(3,1,1,1),mar=c(1,4,0,1))
grow.lf<-Gr(xb=seq(1,k.lf,len=200),r=r.lf,K=k.lf)
harv.lf<-sapply(seq(1,k.lf,len=200),harv,a=a.lf,b=b.lf,c=c.lf,q=q.lf)
plot(y=grow.lf,x=seq(0,1,length.out = length(grow.lf)),ylab="Growth, harvest",type="l",lty=1,lwd=2,cex.lab=1,col="green",xaxt="n",xlab="")
	points(y=harv.lf,x=seq(0,1,length.out = length(harv.lf)),type="l",lwd=2,col="magenta")
	abline(v=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf)/k.lf,lty=3,col="black",lwd=1)
plot(y=yoko(n=seq(1,k.lf,len=200),u=.5,b=0.1,Y=k.lf,M=m.lf)/m.lf,x=seq(0,1,length.out = length(grow.lf)),lty=1, type="l",lwd=2,main="",cex.lab=1,xlab="Population Biomass",col="saddlebrown",
	ylab="Ecosystem service value",cex=.8)
	abline(v=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf)/k.lf,lwd=1,lty=3,col="black")
	arrows(x1=c(-1,-1,-1),x0=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf)/k.lf,
		y0=yoko(n=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf),u=.5,b=.1,Y=k.lf,M=m.lf)/m.lf,lty=3,length=0,col="black",lwd=1)


##########################################
##########################################
##~~~~~~~~~~Compiled figures~~~~~~~~~~~~~~
##########################################
##########################################


###combined growth and ES~population
dev.new(width=6, height=3)
par(mfcol=c(2,3),oma=c(4,0,1,0),cex.axis=1.4,cex.lab=1.5)
#RC
par(mar=c(.5,5,.5,1))
plot(y=grow.cf/k.cf,x=seq(0,1,length.out = length(grow.cf)),ylab="Growth, harvest",type="l",lty=1,lwd=3,
	col="green",xaxt="n",xlab="")
	points(y=harv.cf/k.cf,x=seq(0,1,length.out = length(grow.cf)),type="l",lwd=3,col="magenta")
	abline(v=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf)/k.cf,lty=3,col="black",lwd=3)
text("a", x=.05,y=.043,cex=1.5)
plot(y=yoko(n=seq(1,k.cf,len=200),u=0,b=.1,Y=k.cf,M=m.cf)/m.cf,x=seq(0,1,length.out = length(grow.cf)),lty=1, type="l",lwd=3,main="",
	xlab="Population Biomass",col="sandybrown",ylab="ES value",xaxt="n",yaxt="n",ylim=c(0,1.2))
	axis(side=2,at=c(0,.5,1),labels=c(0,0.5,1))
	axis(side=2,at=c(.25,.75),labels=F)	
	abline(v=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf)/k.cf,lwd=3,lty=3,col="black")
	arrows(x1=c(-1,-1,-1),x0=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf)/k.cf,lwd=3,
		y0=yoko(n=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf),u=0,b=.1,Y=k.cf,M=m.cf)/m.cf,lty=3,length=0,col="black")
	axis(side=1,at=c(0,.5,1))
	axis(side=1,at=c(.25,.75),labels=F)
text("d", x=.05,y=1.15,cex=1.5)
##AC	
par(mar=c(.5,3,.5,3))
plot(y=grow.ac/k.ac,x=seq(0,1,length.out = length(grow.ac)),ylab="",type="l",lty=1,lwd=3,cex.lab=1,col="green",
		xaxt="n",xlab="")
	points(y=harv.ac/k.ac,x=seq(0,1,length.out = length(harv.ac)),type="l",lwd=3,col="magenta")
	abline(v=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac)/k.ac,lty=3,col="black",lwd=3)
text("b", x=.05,y=.17,cex=1.5)
plot(y=yoko(n=seq(1,k.ac,length.out=200),u=1,b=0.1,Y=k.ac,M=m.ac)/m.ac,x=seq(0,1,length.out = length(grow.ac)),lty=1, type="l",lwd=3,main="",
	xlab="",col="sandybrown",ylab="",xaxt="n",yaxt="n",ylim=c(0,1.2))
	abline(v=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac)/k.ac,lwd=3,lty=3,col="black")
	arrows(x1=c(-1,-1,-1),x0=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac)/k.ac,
		y0=yoko(n=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac),u=1,b=.1,Y=k.ac,M=m.ac)/m.ac,lty=3,length=0,col="black",lwd=3)
	axis(side=1,at=c(0,.5,1))
	axis(side=1,at=c(.25,.75),labels=F)
text("e", x=.05,y=1.15,cex=1.5)
#LF
par(mar=c(.5,1,.5,5))
plot(y=grow.lf/k.lf,x=seq(0,1,length.out = length(grow.lf)),ylab="",type="l",lty=1,lwd=3,
		col="green",xaxt="n",xlab="",yaxt="n")
	axis(side=2, at=c(0.005,0.025),labels=c(0.005,0.025))
	points(y=harv.lf/k.lf,x=seq(0,1,length.out = length(harv.lf)),type="l",lwd=3,col="magenta")
	abline(v=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf)/k.lf,lty=3,col="black",lwd=3)
text("c", x=.1,y=.026,cex=1.5)
plot(y=yoko(n=seq(1,k.lf,len=200),u=.5,b=0.1,Y=k.lf,M=m.lf)/m.lf,x=seq(0,1,length.out = length(grow.lf)),lty=1, type="l",lwd=3,main="",
	xlab="Population Biomass",col="sandybrown",ylab="",xaxt="n",yaxt="n",ylim=c(0,1.2))
	abline(v=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf)/k.lf,lwd=3,lty=3,col="black")
	arrows(x1=c(-1,-1,-1),x0=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf)/k.lf,
		y0=yoko(n=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf),u=.5,b=.1,Y=k.lf,M=m.lf)/m.lf,lty=3,length=0,col="black",lwd=3)
	axis(side=1,at=c(0,.5,1))
	axis(side=1,at=c(.25,.75),labels=F)
text("f", x=.05,y=1.15,cex=1.5)
mtext("% Carrying capacity",side=1,cex.lab=1.5,line=2.75,outer=T)


##################
##Population response to Open access
##################

cf.Z<-apply(expand.grid(aa=c(a.cf,a.cf/2,a.cf*2),bb=seq(0,k.cf/5,length.out=200)),1,function (X) 
	c(X[1],X[2]*b.cf,z=Poly(a=X[1],b=X[2],c=c.cf,q=q.cf,r=r.cf,K=k.cf)))
	cf.open.mat<-matrix(data = cf.Z[3,], nrow = 3, ncol = 200)
	matplot(t(cf.open.mat))
ac.Z<-apply(expand.grid(aa=c(a.ac,a.ac*2,a.ac*10),bb=seq(0,k.ac*10,length.out=200)),1,function (X) 
	c(X[1],X[2],z=Poly(a=X[1],b=X[2],c=c.ac,q=q.ac,r=r.ac,K=k.ac)))
	ac.open.mat<-matrix(data = ac.Z[3,], nrow = 3, ncol = 200)
	matplot(t(ac.open.mat))
lf.Z<-apply(expand.grid(aa=c(a.lf,a.lf/5,a.lf*2),bb=seq(0,k.lf/10,length.out=200)),1,function (X) 
	c(X[1],X[2],z=Poly(a=X[1],b=X[2],c=c.lf,q=q.lf,r=r.lf,K=k.lf)))
	lf.open.mat<-matrix(data = lf.Z[3,], nrow = 3, ncol = 200)
	matplot(t(lf.open.mat))

dev.new(width=6, height=3.6)
par(mfcol=c(2,3),oma=c(4,0,0,0),cex.axis=1.4,cex.lab=1.5)
###CF
par(mar=c(1.5,5,1,1))
matplot(y=t(cf.open.mat)/k.cf,x=seq(0,k.cf/5,length.out=200)/k.cf,type="l",lty=c(1,2,3),ylab="% Carrying capacity",xlab="",
	xaxt="n",yaxt="n",col="green",lwd=2,ylim=c(0,1.1))
	axis(side=2,at=c(0,.5,1),labels=c(0,0.5,1))
	axis(side=2,at=c(.25,.75),labels=F)
	points(y=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf)/k.cf,x=b.cf/k.cf,pch=23,cex=2,bg="magenta")
	text("a",x=0.01,y=1.09,cex=1.5,col="black")
matplot(y=t(yoko(n=cf.open.mat,u=0,b=.1,Y=k.cf,M=m.cf)/m.cf),x=seq(0,k.cf/5,length.out=200)/k.cf,yaxt="n",
	type="l",col="sandybrown",lty=1:3,lwd=2,ylab="% Recovery",xlab="",xaxt="s",ylim=c(0,1.1))
	axis(side=2,at=c(0,.5,1),labels=c(0,0.5,1))
	axis(side=2,at=c(.25,.75),labels=F)
	#axis(side=1,at=c(0,.5,1))
	#axis(side=1,at=c(.25,.75),labels=F)
	points(y=yoko(n=Poly(a=a.cf,b=b.cf,c=c.cf,q=q.cf,r=r.cf,K=k.cf),u=0,b=.1,Y=k.cf,M=m.cf)/m.cf 
		,x=b.cf/k.cf,pch=23,cex=2,bg="magenta")
	text("d",x=0.01,y=1.09,cex=1.5,col="black")
###AC
par(mar=c(1.5,3,1,3))
ac.open.mat[ac.open.mat==-Inf]<-NA
matplot(y=t(ac.open.mat)/k.ac,x=seq(0,k.ac*10,length.out=200)/k.ac,type="l",col="green",lty=1:3,lwd=2,ylab="",xlab="",xaxt="n",
	ylim=c(0,1.1),yaxt="n")
	points(y=Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac)/k.ac,x=b.ac/k.ac,pch=23,cex=2,bg="magenta")
	text("b",x=.5,y=1.09,cex=1.5,col="black")
matplot(t(yoko(n=  ac.open.mat,u=1,b=0.1,Y=k.ac,M=m.ac)/m.ac),x=seq(0,k.ac*10,length.out=200)/k.ac,ylim=c(0,1.1),
		type="l",col="sandybrown",lty=1:3,lwd=2,ylab="",xlab="",xaxt="n",yaxt="n")
	axis(side=1,at=c(0,5,10))
	axis(side=1,at=c(1:4,6:9),labels=F)
	points(y=yoko(n= Poly(a=a.ac,b=b.ac,c=c.ac,q=q.ac,r=r.ac,K=k.ac),u=1,b=.1,Y=k.ac,M=m.ac)/m.ac,
		x=b.ac/k.ac,pch=23,cex=2,bg="magenta")
	text("e",x=.5,y=1.09,cex=1.5,col="black")
###LF
par(mar=c(1.5,2,1,4))
matplot(y=t(lf.open.mat)/k.lf,type="l",x=seq(0,k.lf/10,length.out=200)/k.lf,col="green",lty=1:3,lwd=2,
		xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1.1))
	points(y=Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf)/k.lf,x=b.lf/k.lf,pch=23,cex=2,bg="magenta")
	text("c",x=.005,y=1.09,cex=1.5,col="black")
matplot(y=t(yoko(n=lf.open.mat,u=.5,b=0.1,Y=k.lf,M=m.lf))/m.lf,x=seq(0,k.lf/10,length.out=200)/k.lf,
		type="l",col="sandybrown",lty=1:3,lwd=2,xlab="",ylab="",yaxt="n",xaxt="s",ylim=c(0,1.1))
	#axis(side=1,at=c(0,.5,1))
	#axis(side=1,at=c(.25,.75),labels=F)
	points(y=yoko(n= max(Poly(a=a.lf,b=b.lf,c=c.lf,q=q.lf,r=r.lf,K=k.lf),na.rm=T),u=.5,b=.1,Y=k.lf,M=m.lf)/m.lf,
		x=b.lf/k.lf,pch=23,cex=2,bg="magenta")
	text("f",x=.005,y=1.09,cex=1.3,col="black")
mtext("Quantity demanded (by % Carrying capacity)",side=1,line=2,cex.lab=1.3,outer=T)



##################
##Amenity, we picked illistrative values of m
##################

cf.am.Z<-apply(expand.grid(mm=c(m.cf,m.cf*.171,m.cf*0.0000001),bb=seq(1,k.cf/5,length.out=200)),1,function (X)
	c(X[1],X[2],z=social.opt(a=a.cf,r=r.cf,K=k.cf,b=X[2],c=c.cf,q=q.cf,d=d1,y=1,m=X[1])$maximum))
cf.amenity.mat<-matrix(data = cf.am.Z[3,], nrow = 3, ncol = 200)
matplot(t(cf.amenity.mat)/k.cf,ylim=c(0,1))

ac.am.Z<-apply(expand.grid(mm=c(m.ac,m.ac*.01,m.ac*0.0000001),bb=seq(1,k.ac*10,length.out=200)),1,function (X)
	c(X[1],X[2],z=social.opt(a=a.ac,r=r.ac,K=k.ac,b=X[2],c=c.ac,q=q.ac,d=d1,y=4,m=X[1])$maximum))
ac.amenity.mat<-matrix(data = ac.am.Z[3,], nrow = 3, ncol = 200)
matplot(t(ac.amenity.mat)/k.ac,ylim=c(0,1))

lf.am.Z<-apply(expand.grid(mm=c(m.lf,m.lf*.05,m.lf*0.0000001),bb=seq(1,k.lf/10,length.out=200)),1,function (X)
	c(X[1],X[2],z=social.opt(a=a.lf,r=r.lf,K=k.lf,b=X[2],c=c.lf,q=q.lf,d=d1,y=2,m=X[1])$maximum))
lf.amenity.mat<-matrix(data = lf.am.Z[3,], nrow = 3, ncol = 200)
matplot(t(lf.amenity.mat)/k.lf,ylim=c(0,1))


dev.new(width=6, height=3.6)
par(mfcol=c(2,3),oma=c(4,0,0,0),cex.axis=1.4,cex.lab=1.5)
###CF
par(mar=c(1.5,5,1,1))
matplot(y=t(cf.amenity.mat)/k.cf,x=seq(0,k.cf/5,length.out=200)/k.cf,type="l",lty=c(1,2,3),ylab="% Carrying capacity",xlab="",
	xaxt="n",yaxt="n",col="green",lwd=2,ylim=c(0,1.1))
	axis(side=2,at=c(0,.5,1),labels=c(0,0.5,1))
	axis(side=2,at=c(.25,.75),labels=F)
	text("a",x=.01,y=1.09,cex=1.5,col="black")
	#points(y= social.opt(a=a.cf,r=r.cf,K=k.cf,b=b.cf,c=c.cf,q=q.cf,d=d1,y=1,m=m.cf)$maximum/k.cf ,x=b.cf/k.cf,  pch=23,cex=2,bg="magenta")
matplot(y=t(yoko(n=cf.amenity.mat,u=0,b=.1,Y=k.cf,M=m.cf)/m.cf),x=seq(0,k.cf/5,length.out=200)/k.cf,xaxt="s",yaxt="n",
	type="l",col="sandybrown",lty=1:3,lwd=2,ylab="% Recovery",xlab="",ylim=c(0,1.1))
	axis(side=2,at=c(0,.5,1),labels=c(0,0.5,1))
	axis(side=2,at=c(.25,.75),labels=F)
	#axis(side=1,at=c(0,.5,1),labels=c(0,0.5,1))
	#axis(side=1,at=c(.25,.75),labels=F)
	text("d",x=0.01,y=1.09,cex=1.5,col="black")
	#points(y= yoko(n=social.opt(a=a.cf,r=r.cf,K=k.cf,b=b.cf,c=c.cf,q=q.cf,d=d1,y=1,m=m.cf)$maximum,u=0,b=.1,Y=k.cf,M=m.cf)/m.cf,
	#x=b.cf/k.cf,  pch=23,cex=2,bg="magenta")
###AC
par(mar=c(1.5,3,1,3))
ac.amenity.mat[ac.amenity.mat==-Inf]<-NA
matplot(y=t(ac.amenity.mat)/k.ac,x=seq(0,k.ac*10,length.out=200)/k.ac,type="l",col="green",lty=1:3,lwd=2,ylab="",xlab="",xaxt="n",
	ylim=c(0,1.1),yaxt="n")
	text("b",x=1,y=1.09,cex=1.5,col="black")
	#points(y= social.opt(a=a.ac,r=r.ac,K=k.ac,b=b.ac,c=c.ac,q=q.ac,d=d1,y=1,m=m.ac)$maximum/k.ac ,x=b.ac/k.ac,  pch=23,cex=2,bg="magenta")
matplot(t(yoko(n=  ac.amenity.mat,u=1,b=0.1,Y=k.ac,M=m.ac)/m.ac),x=seq(0,k.ac*10,length.out=200)/k.ac,yaxt="n",
		type="l",col="sandybrown",lty=1:3,lwd=2,ylab="",xlab="",ylim=c(0,1.1),xaxt="n")
	axis(side=1,at=c(0,5,10))
	axis(side=1,at=c(1:4,6:9),labels=F)
	text("e",x=.5,y=1.09,cex=1.5,col="black")
	#points(y= yoko(n=social.opt(a=a.ac,r=r.ac,K=k.ac,b=b.ac,c=c.ac,q=q.ac,d=d1,y=1,m=m.ac)$maximum,u=1,b=.1,Y=k.ac,M=m.ac)/m.ac,
	#x=b.ac/k.ac,  pch=23,cex=2,bg="magenta")

###LF
par(mar=c(1.5,2,1,4))
matplot(y=t(lf.amenity.mat)/k.lf,type="l",x=seq(0,k.lf/10,length.out=200)/k.lf,col="green",lty=1:3,lwd=2,xlab="",
	ylab="",xaxt="n",yaxt="n",,ylim=c(0,1.1))
	text("c",x=.01,y=1.09,cex=1.5,col="black")
	#points(y= social.opt(a=a.lf,r=r.lf,K=k.lf,b=b.lf,c=c.lf,q=q.lf,d=d1,y=1,m=m.lf)$maximum/k.lf ,x=b.lf/k.lf,  pch=23,cex=2,bg="magenta")
matplot(y=t(yoko(n=lf.amenity.mat,u=.5,b=0.1,Y=k.lf,M=m.lf))/m.lf,x=seq(0,k.lf/10,length.out=200)/k.lf,xaxt="s",
	type="l",col="sandybrown",lty=1:3,lwd=2,xlab="",ylab="",yaxt="n",ylim=c(0,1.1))
	#axis(side=1,at=c(0,.5,1),labels=c(0,0.5,1))
	#axis(side=1,at=c(.25,.75),labels=F)
	text("f",x=.01,y=1.09,cex=1.3,col="black")
	#points(y= yoko(n=social.opt(a=a.lf,r=r.lf,K=k.lf,b=b.lf,c=c.lf,q=q.lf,d=d1,y=1,m=m.lf)$maximum,u=1,b=.1,Y=k.lf,M=m.lf)/m.lf,
	#x=b.lf/k.lf,  pch=23,cex=2,bg="magenta")
mtext("Quantity demanded (by % Carrying capacity)",side=1,line=2,cex.lab=1.3,outer=T)

######################
#####~~~~~~~~subsidy
######################

cf.am.sub<-apply(expand.grid(mm=c(m.cf,m.cf*.171,m.cf*0.0000001),bb=seq(1,k.cf/5,length.out=200)),1,function (X)
	c(X[1],X[2],z=soc.sub(a=a.cf,r=r.cf,K=k.cf,b=X[2],c=c.cf,q=q.cf,y=1,m=X[1])))
cf.sub.mat<-matrix(data = cf.am.sub[3,], nrow = 3, ncol = 200)
matplot(sapply(1:3,function (X) cf.sub.mat[X,]/c(m.cf,m.cf*.07,m.cf*0.0000001)[X])-1)

ac.am.sub<-apply(expand.grid(mm=c(m.ac,m.ac*.01,m.ac*0.0000001),bb=seq(1,k.ac*10,length.out=200)),1,function (X)
	c(X[1],X[2],z=soc.sub(a=a.ac,r=r.ac,K=k.ac,b=X[2],c=c.ac,q=q.ac,y=4,m=X[1])))
ac.sub.mat<-matrix(data = ac.am.sub[3,], nrow = 3, ncol = 200)
matplot(sapply(1:3,function (X) ac.sub.mat[X,]/c(m.ac,m.ac*.01,m.ac*0.0000001)[X])-1)

lf.am.sub<-apply(expand.grid(mm=c(m.lf,m.lf*.05,m.lf*0.0000001),bb=seq(1,k.lf/10,length.out=200)),1,function (X)
	c(X[1],X[2],z=soc.sub(a=a.lf,r=r.lf,K=k.lf,b=X[2],c=c.lf,q=q.lf,y=2,m=X[1])))
lf.sub.mat<-matrix(data = lf.am.sub[3,], nrow = 3, ncol = 200)
matplot(sapply(1:3,function (X) lf.sub.mat[X,]/c(m.lf,m.lf*.05,m.lf*0.0000001)[X])-1)


dev.new(width=6, height=2.1)
par(mfcol=c(1,3),oma=c(4,1,0,0),cex.axis=1.4,cex.lab=1.5)
###CF
par(mar=c(1.5,5,1,1))
matplot(y=sapply(1:3,function (X) cf.sub.mat[X,]/c(m.cf,m.cf*.171,m.cf*0.0000001)[X])-1,x=seq(0,k.cf/5,length.out=200)/k.cf,type="l",
	lty=c(1,2,3),ylab="% ES value",xlab="",xaxt="s",yaxt="n",col="magenta",lwd=2,ylim=c(0,250000))
	axis(side=2,at=c(0,100000,200000),labels=c("0","1e5","2e5"))
	axis(side=2,at=c(50000,150000,250000),labels=F)
	#axis(side=1,at=c(0,.5,1))
	#axis(side=1,at=c(.25,.75),labels=F)
	text("a",x=.01,y=240000,cex=1.5,col="black")
###AC
par(mar=c(1.5,4,1,2))
ac.amenity.mat[ac.amenity.mat==-Inf]<-NA
matplot(y=sapply(1:3,function (X) ac.sub.mat[X,]/c(m.ac,m.ac*.01,m.ac*0.0000001)[X])-1,x=seq(0,k.ac*10,length.out=200)/k.ac,
	type="l",col="magenta",lty=1:3,lwd=2,ylab="",xlab="",xaxt="n",yaxt="s",ylim=c(-90,10))
	axis(side=1,at=c(0,5,10))
	axis(side=1,at=c(1:4,6:9),labels=F)
	text("b",x=1,y=7,cex=1.5,col="black")
###LF
par(mar=c(1.5,2,1,4))
matplot(y=sapply(1:3,function (X) lf.sub.mat[X,]/c(m.lf,m.lf*.05,m.lf*0.0000001)[X])-1,type="l",x=seq(0,k.lf/15,length.out=200)/k.lf,
	col="magenta",lty=1:3,lwd=2,xlab="",ylab="",xaxt="s",yaxt="n",ylim=c(0,250000))
	#axis(side=1,at=c(0,.5,1))
	#axis(side=1,at=c(.25,.75),labels=F)
	axis(side=2,at=c(0,100000,200000),labels=c("0","1e5","2e5"))
	axis(side=2,at=c(50000,150000,250000),labels=F)
	text("c",x=.005,y=240000,cex=1.5,col="black")
mtext("Quantity demanded (by % Carrying capacity)",side=1,line=2,cex.lab=1.3,outer=T)
mtext("Social welfare as",side=2,line=-.50,cex.lab=1.3,outer=T)






#########################################