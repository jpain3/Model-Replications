#Model Replication
#Adaptation as Information Restriction
#Jerker Denrell, James G. March, 2001

#Created by James Paine
#15.879 - Simulation Models in Social and Behavioral Sciences
#February 22, 2019

#####PART 2 - EXPERIENTIAL LEARNING WITH COMPETENCY BUILDING#####


###System Parameters

#Expectation of the 'risky' or 'new' alternative
X = 15
#Standard deviation of the 'risky' or 'new' alternative
S = 5

#Expectation of the 'certain' alternative
Y = 10

#Rate at which aspirations adjust to experience
b= 0.5

#Initial Competence of agents with the 'risky alternative
c0 = 0.3

#Rate of standard deviation decrease as a function of competence (see page 527)
k = 0.5


###Simulation Paramters
periods = 50
agents = 5000
replications = 5

#Display full iteration details - adds time, only use for small replications
fulldetail = 0

#Create arrays to range over for later plotting
A = seq(from = .05, to = 0.95, by = ((0.95-0.05)/18))
D = c(0.1,0.3)

#Note P is the probability that the RISKY alternative is chosen
#This is a two-choice system, so the probability of the certain
#  choice is 1-P for any timestep t


#Note on array indicides
#P[n,t,a] and L[n,t,a] and choice[n,t,a] where
# n = agent number
# t = period number
# P = probability of choosing the risky choice
# example: P[50,20] = 0.25 means that agent #50 had a 25% chance of 
#   choosing the risky choice at time period 20

#IN the below sampling, choice 1 is the 'Risky' alternative

#Initialize the probability array P and the aspiration array L
#initalize the GLOBAL probability and q array spaces
  P = array(data = NA, dim = c(agents,periods+1,length(A),length(D),replications))
  L = array(data = NA, dim = c(agents,periods+1,length(A), length(D),replications))
  choice = array(data = NA, dim = c(agents,periods+1,length(A), length(D),replications))
  FractX = array(data = NA, dim = c(periods+1,length(A), length(D),replications))
  c = array(data = NA, dim = c(agents,periods+1,replications))
##Initialize the simluation

  #Initialize the probabilities at t=0
  P[,1,,,]=0.5
  
  #Initialize the asipriations at t=0
  L[,1,,,] = P[,1,,,]*X+(1-P[,1,,,])*Y

  #Initialize the competencies at t=0
  c[,1,] = c0
  
n=1
t=1
SoL = 1
LearnParm = 1
SoL = 1

starttime = proc.time()

for (rep in 1:replications) {
  
  if (fulldetail != 1) {
    cat("\014")
    print(paste0("Replication ",rep," of ", replications))
    print(paste0("Elapsed time since last replication: ",(proc.time() - starttime)[3],"s"))
    print(paste0("Avg Replication time: ",(proc.time() - starttime)[3]/rep,"s per replication"))
    print(paste0("Est Time Remaining: ",((proc.time() - starttime)[3]/rep)*(replications-rep),"s"))
    
    
  }

  for (LearnParm in 1:length(D)) {
    
    d = D[LearnParm]
  
    for (SoL in 1:length(A)) {  #Note: a is the 'speed of learning'
    
      a = A[SoL]
      
      for (t in 1:(periods+1)) {
        
        if (fulldetail == 1) {
          cat("\014")
          print(paste0("Replication ",rep," of ", replications))
          print(paste0("d value: ", d))
          print(paste0("a value: ", a))
          print(paste0("Time Period: ", (t-1), " of ", periods))
        }
        
        for (n in 1:agents) {
          
          #print(paste0("Agent: ", n, " of ", agents))
          
          #Probabilistically get agent n's choice
          # Note, here choice 1 is the 'risky' or 'new' alternative
          choice[n,t,SoL,LearnParm,rep]=sample(c(1:2),1,replace=FALSE,prob=c(P[n,t,SoL,LearnParm,rep],(1-P[n,t,SoL,LearnParm,rep])))
          
          #Get the realized value from either the risky or certain distributions
          if (choice[n,t,SoL,LearnParm,rep] == 1) {
            
            #Get the st dev based on the current competency c
            stdev = (S/c[n,t,rep])^k
            avg = c[n,t,rep]*X
            #Get the performance experienced
            O = rnorm(1, mean = avg, sd = stdev)
          } else {
            O = Y
          }
          
          
          if ((t+1)<=(periods+1)) {
          
            #Get the aspirations for next round
            L[n,t+1,SoL,LearnParm,rep] = L[n,t,SoL,LearnParm,rep]*(1-b)+O*b
            
            #Update the probability of choosing the risky alternative based on the
            #   current aspiration L and the realized value O
            
            P[n,t+1,SoL,LearnParm,rep] = P[n,t,SoL,LearnParm,rep]
            
            #"...if the first/new alternative is tried..."
            if (choice[n,t,SoL,LearnParm,rep] == 1) {
              
              #"...Competence increases with each utilization"
              c[n,t+1,rep]=c[n,t,rep] + d*(1-c[n,t,rep])
  
              #"...AND yields an outcome better than the asipration at time t"
              if (O > L[n,t,SoL,LearnParm,rep]) {
                #"The probability increases in the following way:"
                P[n,t+1,SoL,LearnParm,rep]=P[n,t,SoL,LearnParm,rep] + a*(1-P[n,t,SoL,LearnParm,rep])
                
              #"...OR yields an outcome that is worse than the asipriation
              } else if (O < L[n,t,SoL,LearnParm,rep]) {
                P[n,t+1,SoL,LearnParm,rep] = (1-a)*P[n,t,SoL,LearnParm,rep]
              
              }
            }
            
            #...or if the second/existing alternative is tried..."
            if (choice[n,t,SoL,LearnParm,rep] == 2) {
              
              c[n,t+1,rep]=c[n,t,rep]
              
              #"...AND yields an outcome worse than the asipration at time t"
              if (O < L[n,t,SoL,LearnParm,rep]) {
                #"The probability increases in the following way:"
                P[n,t+1,SoL,LearnParm,rep]=P[n,t,SoL,LearnParm,rep] + a*(1-P[n,t,SoL,LearnParm,rep])
                
                #"...OR yields an outcome that is better than the asipriation
              } else if (O > L[n,t,SoL,LearnParm,rep]) {
                P[n,t+1,SoL,LearnParm,rep] = (1-a)*P[n,t,SoL,LearnParm,rep]
                
              }        
              
            }
      
          }
          
          
        } # Next Agent n 
    
        # Determine fractions of choices at each time step (t) and Speed of Learning (a) value
        NumX = sum(choice[,t,SoL,LearnParm,rep] == 1)
        NumY = sum(choice[,t,SoL,LearnParm,rep] == 2)
        FractX[t,SoL,LearnParm,rep] = NumX/(NumX+NumY)
        
            
      } # Next time period t
    
    
    } # Next speed of learning parameter a value
  
  } # Next learning paramter d value

} # Next replication


#Average each replication and determine the standard deviation
FractX_avg = apply(FractX, c(1,2,3), mean)
FractX_sd = apply(FractX, c(1,2,3), sd)

#Plot specific repitition and round combo
rep = 1
Round = 50
plot(A,FractX[Round,,1,rep], typ = "l", col = "black", lwd=2, lty =1,
     xlim=c(0,1), ylim=c(0,1),
     xlab = "Speed of Learning Paramter (a)", ylab = "Fraction of Agents Choosing Risking Alternative")
  points(A,FractX[Round,,1,rep], typ = "p", pch=4, lwd=2)
  lines(A,FractX[Round,,2,rep], type = "l", col = "black", lwd=2, lty=1)
  points(A,FractX[Round,,2,rep], typ = "p", pch=0, lwd=2)
  text(c(0.21,0.3),c(0.25,0.48),paste("d = ", D), cex = .75)
  
  
##FIGURE 2

Round = 50
plot(A,FractX_avg[Round,,1], typ = "l", col = "black", lwd=2, lty =3,
     xlim=c(0,1), ylim=c(0,1),
     xlab = "Speed of Learning Paramter (a)", ylab = "Fraction of Agents Choosing Risking Alternative")
points(A,FractX_avg[Round,,1], typ = "p", pch=4, lwd=2)
lines(A,FractX_avg[Round,,2], type = "l", col = "black", lwd=2, lty=1)
points(A,FractX_avg[Round,,2], typ = "p", pch=22, lwd=2, bg = "white")
text(c(0.21,0.3),c(0.25,0.48),paste("d = ", D), cex = .75)

