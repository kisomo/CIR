#Calibrating and Projecting the Cox-Ingersoll-Ross Process
#ACSC 351 Final Project
#Matt Galloway and Caleb Moore

#note: a "#" designates a line of code as a comment
#Before running this script, you must specify your 'working directory'. This is the location where any output will be stored.
#For Mac: 'MISC'->'Change Working Directory...' (And then specify location)
#For PC: Exectue the following line of code: setwd('location'). Excute by "CNTRL+R". Example below...
# setwd("C:/Users/gall3017/Desktop/CIR")

#The only thing you must specify (besides working directory) is the 'ncap' - or number of periods you wish to forecast - and the 'sim' - the number of simulations you wish to run

#To run file...
#For Mac: Highlight all of the lines of code and then 'Edit'->'Execute'
#For PC: 'Edit'->'Run All'
#note: You will then be prompted to choose your data file in the browser. Once you specify the file, R will then continue to execute the remaining lines of code. Your data file must be saved as a ".csv" and interest rates you wish to calibrate/forecast must be named 'Data' in the column header.

#Your ouput will be saved as a ".png" and ".txt" file in your working directory.

#Equation...
#CIR: dr=a*(b-r)*dt+sigma*sqrt(r)dZ


#delete previously stored values
rm(list=ls())

setwd("C:/Users/gall3017/Desktop")


#Be Consistent! (If you are using WEEKLY data, make sure your ncap is in WEEKS)

#number of simulations
sim<-50
#Single premium payment
premium<-100000
#Length of annuity
years<-15
#Payments per year
paymentsper<-12


#Total payments
ncap<-years*paymentsper



#choose data file
InterestRates<-read.csv(file.choose())
#display structure of data
str(InterestRates)
#summary of data
summary(InterestRates)
#attach names
attach(InterestRates)
#display names
names(InterestRates)
#download 'lattice' package from library
library(lattice)


#CALIBRATION


#data
Data<-InterestRates$Data
#plot the data
#plot(Counter,InterestRates$X5.year, xlab="Week", ylab="Interest Rate (5yr.)", type='l',col='red')
#title(main="5yr. Interest Rates (1/1/02-11/20/13)", col.main="black", font.main=2)

#Da - discrete drift term
#Dsigma - discrete volatility term
#a - continuous drift term
#b - long-run average term
#sigma - continuous volatility term
#N - number of data values


#calculate number of entries
N<-length(InterestRates$Data)
#calculate long-run average
b<-mean(InterestRates$Data)


#designate function that will be used to calibrate the parameters
#more specifically, this is the sum of residuals (RSS)
calibration<-function(param){sum(((Data[2:N]-b)-param*(Data[1:N-1]-b))^2/((Data[2:N]-b)+b))}
#minimize (nlm) the residuals
results<-nlm(calibration,param<-0.9901,hessian=TRUE) #minimization
#display results
results


#calculate discrete drift
Da<-results$estimate
#calculate discrete volatility
Dsigma<-sqrt(results$minimum/(N-1))
#calculate continuous drift, given discrete drift
a<-log(1/results$estimate)
#calculate continuous volatility, given discrete volatility
sigma<-sqrt((Dsigma^2)*2*a/(1-exp(-2*a)))
#calculate minimum RSS
RSS<-results$minimum
#expected value of forecasted interest rate
E<-(b+(InterestRates$Data[N]-b)*exp(-a*ncap))
#variance of forecasted interest rate
V<-(((InterestRates$Data[N]*sigma^2*exp(-a*ncap))/a)*(1-exp(-a*ncap))+(b*sigma^2/2/a)*(1-exp(-a*ncap))^2)


#SIMULATION
 
SIMULATION=function(t.start,x.start){
	X<-array(0,c(ncap,sim))
	for(j in 1:sim){
	t.old=t.start
	x.old=x.start
	d=4*a*b/sigma^2
	for(i in 1:ncap){
		t.new=t.old+1
		n=(4*a*exp(-a*(t.new-t.old)))/(sigma^2*(1-exp(-a*(t.new-t.old))))
		P=rpois(1,0.5*x.old*n)
		v=d+2*P
		Z=rchisq(1,v,ncp=0)
		x.new=Z*exp(-a*(t.new-t.old))/n
		X[i,j]=x.new
		x.old=x.new
	}
	}
return(list(X=X))
}

#set seed for random generation
set.seed(400)
#set initial values
simulate=SIMULATION(t.start=0,x.start=InterestRates$Data[N])


#designate a new 'Counter' variable to account for added projection values
new.Counter<-c(1:length(c(InterestRates$Data[1:N],simulate$X[,1])))


#PRICING

guaranteed_interest<-((1+InterestRates$Data[N]/100)^(1/paymentsper)-1)
payment<-((premium*guaranteed_interest)/(1-(1/(1+guaranteed_interest))^(ncap)))

PRICING=function(facevalue.start){
	M<-array(0,c(ncap,sim))
	for(k in 1:sim){
	facevalue.old=facevalue.start
	for(l in 1:ncap){
		facevalue.new=facevalue.old*((1+simulate$X[l,k]/100)^(1/paymentsper))-payment
		M[l,k]=facevalue.new
		facevalue.old=facevalue.new
		}
	}
return(list(M=M))
}

#set initial values
price=PRICING(facevalue.start=premium)


#GRAPHING

#'print' plot output to a .png file located in workspace
#designate file name, size, resolution
png(file="CIR_Pricing_Projections(monthly).png",width=11,height=8,units="in",res=600)
par(mfrow=c(2,1))
for(j in 1:sim){
#plot projections, name axis, color
plot(new.Counter,c(InterestRates$Data[1:N],simulate$X[,j]), xlab="", ylab="Interest Rate",lwd=1,ylim=c(0,max(c(InterestRates$Data,simulate$X))+0.1),type='l',col='black');
par(new=TRUE)}
#include vertical line showing where original data ends and projections begin
lines(new.Counter,c(InterestRates$Data[1:N],rowMeans(simulate$X,na.rm=FALSE,dims=1)),type='l',lwd=1,col='orange')
abline(v=N,untf=FALSE,col='red',lwd=2)
#add title to graph
title(main="Interest Rates Projections", col.main="black", font.main=2)

for(j in 1:sim){
#plot projections, name axis, color
plot(new.Counter,c(InterestRates$Data[1:N],simulate$X[,j]), xlab="Periods (Months)", ylab="Interest Rate",lwd=1,ylim=c(min(simulate$X)-0.1,max(simulate$X)+0.1),xlim=c(N-0.05*N,length(c(InterestRates$Data[1:N],simulate$X[,1]))),type='l',col='black');
par(new=TRUE)}
#include vertical line showing where original data ends and projections begin
lines(new.Counter,c(InterestRates$Data[1:N],rowMeans(simulate$X,na.rm=FALSE,dims=1)),type='l',lwd=2,col="orange")
abline(v=N,untf=FALSE,col='red',lwd=3)
dev.off()

#"sink" the output to an external data file
sink("CIR_Pricing_Projections(monthly).txt")
print("Cox-Ingersoll-Ross Project")
cat("Discrete Drift Term (Da)",sep="/n")
results$estimate
cat("Discrete Volatility (Dsigma)",sep="/n")
Dsigma
cat("Continuous Drift (a)",sep="/n")
a
cat("Continuous Volatility (sigma)",sep="/n")
sigma
cat("Long-Run Average (b)",sep="/n")
b
cat("RSS")
results$minimum
cat("Number of Entries (N)",sep="/n")
N
print("Guaranteed Rate (Current - Annual)")
InterestRates$Data[N]
print("Single Premium (in dollars)")
premium
print("Level Monthly Payments")
payment
cat("Length of Projection (ncap)(Months)",sep="/n")
ncap
cat("Profit for Insurance Company",sep="/n")
for(k in 1:sim){
	print(price$M[ncap,k])
	}
cat("Average Profit",sep="/n")
rowMeans(price$M,na.rm=FALSE,dims=1)[ncap]