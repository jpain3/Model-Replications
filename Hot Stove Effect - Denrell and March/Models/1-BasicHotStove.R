#Model Replication
#Adaptation as Information Restriction
#Jerker Denrell, James G. March, 2001

#Created by James Paine
#15.879 - Simulation Models in Social and Behavioral Sciences
#February 22, 2019

#####PART 1 - EXPERIENTIAL LEARNING WITHOUT COMPETENCY BUILDING#####


###System Parameters

#Expectation of the 'risky' alternative
X = 10
#Standard deviation of the 'risky alternative
S = 10

#Expectation of the 'certain' alternative
Y = 10

#Rate at which aspirations adjust to experience
b= 0.5


###Simulation Paramters
periods = 50
agents = 5000
replications = 2

A = seq(from = 0, to = 0.045, by = ((0.05-0)/5))
A = append(A,seq(from = 0.05, to = 1, by = ((1-0)/20)))
A

#Note P is the probability that the RISKY alternative is chosen
#This is a two-choice system, so the probability of the certain
#  choice is 1-P for any timestep t


#Note on array indicides
#P[n,t,a] and L[n,t,a] and k[n,t,a] where
# n = agent number
# t = period number
# P = probability of choosing the risky choice
# example: P[50,20] = 0.25 means that agent #50 had a 25% chance of 
#   choosing the risky choice at time period 20

#IN the below sampling, choice 1 is the 'Risky' alternative

#Initialize the probability array P and the aspiration array L
#initalize the GLOBAL probability and q array spaces
  P = array(data = NA, dim = c(agents,periods+1,length(A),replications))
  L = array(data = NA, dim = c(agents,periods+1,length(A),replications))
  k = array(data = NA, dim = c(agents,periods+1,length(A),replications))
  FractX = array(data = NA, dim = c(periods+1,length(A),replications))

##Initialize the simluation

  #Initialize the probabilities at t=0
  P[,1,,]=0.5
  
  #Initialize the asipriations at t=0
  L[,1,,] = P[,1,,]*X+(1-P[,1,,])*Y

#debug param
n=1
t=1
SoL = 1
rep = 1
#

for (rep in 1:replications) {

  for (SoL in 1:length(A)) {
  
    a = A[SoL]
    
    for (t in 1:(periods+1)) {
      
      cat("\014")
      print(paste0("Replication ", rep," of ", replications))
      print(paste0("a value: ", a))
      print(paste0("Time Period: ", (t-1), " of ", periods))
      
      
      for (n in 1:agents) {
        
        #print(paste0("Agent: ", n, " of ", agents))
        
        #Probabilistically get agent n's choice 
        k[n,t,SoL,rep]=sample(c(1:2),1,replace=FALSE,prob=c(P[n,t,SoL,rep],(1-P[n,t,SoL,rep])))
        
        #Get the realized value from either the risky or certain distributions
        if (k[n,t,SoL,rep] == 1) {
          O = rnorm(1, mean = X, sd = S)
        } else {
          O = Y
        }
        
        
        if ((t+1)<=(periods+1)) {
        
          #Get the aspirations for next round
          L[n,t+1,SoL,rep] = L[n,t,SoL,rep]*(1-b)+O*b
          
          #Update the probability of choosing the risky alternative based on the
          #   current aspiration L and the realized value O
          
          P[n,t+1,SoL,rep] = P[n,t,SoL,rep]
          
          #"...if the first alternative is tried..."
          if (k[n,t,SoL,rep] == 1) {
            #"...AND yields an outcome better than the asipration at time t"
            if (O > L[n,t,SoL,rep]) {
              #"The probability increases in the following way:"
              P[n,t+1,SoL,rep]=P[n,t,SoL,rep] + a*(1-P[n,t,SoL,rep])
              
            #"...OR yields an outcome that is worse than the asipriation
            } else if (O < L[n,t,SoL,rep]) {
              P[n,t+1,SoL,rep] = (1-a)*P[n,t,SoL,rep]
            
            }
          }
          
          #...or if the second alternative is tried..."
          if (k[n,t,SoL,rep] == 2) {
            #"...AND yields an outcome worse than the asipration at time t"
            if (O < L[n,t,SoL,rep]) {
              #"The probability increases in the following way:"
              P[n,t+1,SoL,rep]=P[n,t,SoL,rep] + a*(1-P[n,t,SoL,rep])
              
              #"...OR yields an outcome that is better than the asipriation
            } else if (O > L[n,t,SoL,rep]) {
              P[n,t+1,SoL,rep] = (1-a)*P[n,t,SoL,rep]
              
            }        
            
          }
    
        }
        
        
      } # Next Agent n 
  
      # Determine fractions of choices at each time step (t) and Speed of Learning (a) value
      NumX = sum(k[,t,SoL,rep] == 1)
      NumY = sum(k[,t,SoL,rep] == 2)
      FractX[t,SoL,rep] = NumX/(NumX+NumY)
      
          
    } # Next time period t
  
  
  } # Next parameter a value
  
} # Next replication


#Average each replication and determine the standard deviation
FractX_avg = apply(FractX, c(1,2), mean)
FractX_sd = apply(FractX, c(1,2), sd)

#Get Confidence intervals for plotting
CI = 0.999
SD_mod = qnorm(1-CI,lower.tail=FALSE)/sqrt(replications)
FractX_UCI = FractX_avg + SD_mod*FractX_sd
FractX_LCI = FractX_avg - SD_mod*FractX_sd
FractX_UCI[FractX_UCI>1] = 1
FractX_LCI[FractX_LCI<0] = 0


#Full plot of all rounds of last replication
Rep = replications
matplot(A,t(FractX[,,Rep]), typ = "l", pch=4, col = "grey", lwd=1, lty =1,
        xlim=c(0,1), ylim=c(0,1),
        xlab = "Speed of Learning Paramter (a)", ylab = "Fraction of Agents Choosing Risking Alternative"
)


#Plot of last round of specific replication
Round = periods
Rep = 2
plot(A,FractX[Round,,Rep], typ = "l", pch=4, col = "black", lwd=2, lty =1,
     xlim=c(0,1), ylim=c(0,1),
     xlab = "Speed of Learning Paramter (a)", ylab = "Fraction of Agents Choosing Risking Alternative"
)
points(A,FractX[Round,,Rep], typ = "p", pch=4, lwd=2)


#Full plot all rounds of the mean of all replications
matplot(A,t(FractX_avg[,]), typ = "l", pch=4, col = "grey", lwd=1, lty =1,
        xlim=c(0,1), ylim=c(0,1),
        xlab = "Speed of Learning Paramter (a)", ylab = "Fraction of Agents Choosing Risking Alternative"
)

##FIGURE 1 RECREATION
Round = periods
plot(A,FractX_avg[Round,], typ = "l", pch=4, col = "black", lwd=2, lty =1,
     xlim=c(0,1), ylim=c(0,1),
     xlab = "Speed of Learning Paramter (a)", ylab = "Fraction of Agents Choosing Risking Alternative"
)
points(A,FractX_avg[Round,], typ = "p", pch=4, lwd=2)

#arrows(A,FractX_LCI[Round,], A, FractX_UCI[Round,], length=0.05, angle=90, code=3)
