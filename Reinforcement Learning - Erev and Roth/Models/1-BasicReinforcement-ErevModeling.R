#Model Replication
#Predicting How People Play Games: Reinforcement Learning in Experimental Games with Unique, Mixed Strategy Equilibria 
#Erev and Roth, 1998

#Created by James Paine
#15.879 - Simulation Models in Social and Behavioral Sciences
#February 14, 2019


#ONE-PARAMETER MODEL


setwd("C:/Users/jpaine/OneDrive/MIT/15.879 - Simulation Models in Social and Behavioral Sciences/Erev Modeling")
Payout = read.csv("ML Payoff Matrix.csv")
head(Payout)

#Set the number of time steps or rounds of individual play between each player
Rounds = 1000

#Set the number of players
#NOTE: MUST MATCH FORMAT OF INPUT PAYOUT MATRIX
NumPlayer = 2

#Set the number of interations of each play
Iterations = 20

#Set Strength Parameter
#Tuned value on page 862
S1 = 54


#Note on indicies
#i = iteration index
#n = player index
#t = round, or time, index of single play
#k = index of choice the player actually played
#j = general index of choices available to player
#q[n,j,t] = 'propensity' towards a strategy j
#p[n,j,t] = probability of choosing strategy j



#Get the number of unique choices for each player
M = lapply(Payout, function(x) length(unique(x)))
M = M[1:NumPlayer]

#Get the average absolute payoff for each player
X = lapply(Payout, function(x) abs(mean(x)))
X = X[(NumPlayer+1):length(X)]

#Get the minimum payoff for each player
x_min = lapply(Payout, function(x) min(x))
x_min = x_min[(NumPlayer+1):length(x_min)]

#initalize the GLOBAL probability and q array spaces
P = array(data = NA, dim = c(NumPlayer,max(as.numeric(M[1:NumPlayer])),Rounds,Iterations))
Q=P

#initialize the GLOBAL choice array k
K = array(data = NA, dim = c(NumPlayer,Rounds,Iterations))

i=1

#Loop through each iteration i the multi-round game
for (i in 1:Iterations) {
  cat("\014")  
  print(paste0("Iteration: ", i, " of ", Iterations))
  
  
  #initalize the probability and q array space for this iteration
  p = array(data = NA, dim = c(NumPlayer,max(as.numeric(M[1:NumPlayer])),Rounds))
  q=p
  
  #initialize the choice array k for this iteration
  k = array(data = NA, dim = c(NumPlayer,Rounds))
  
  
  #Initialization round
  #set t=1 values for the probabilities and propensities
  for (n in 1:NumPlayer) {
    for (j in 1:as.numeric(M[n])) {
      
      p[n,j,1]= 1/as.numeric(M[n])
      q[n,j,1]= p[n,j,1]*S1*as.numeric(X[n])
    }
  }

  #Loop through subsequent rounds
  for(t in 2:Rounds) {
    
    for(n in 1:NumPlayer){
      #Probabilistically get player n's choice 
      k[n,t]=sample(c(1:as.numeric(M[n])),1,replace=FALSE,prob=p[n,,(t-1)])
    }
    
    #Look at the payout matrix to determine which row matches the choices
    FullPayout = Payout[apply(Payout[,1:NumPlayer], 1, function(x) identical(as.numeric(x),as.numeric(k[,t]))),]
    R = FullPayout[,(NumPlayer+1):ncol(Payout)]-x_min
    
    for(n in 1:NumPlayer){
      
      #for each strategy option, see if a reward was observed
      for(j in 1:as.numeric(M[n])){
        
        q[n,j,t] = q[n,j,(t-1)]
        
        if (as.numeric(k[n,t]) == j) {
          #update the propensities based on the reward observed
          q[n,j,t] =  q[n,j,(t-1)]+as.numeric(R[n])         
        }
      } #next j in the observation of choice and propentsity update loop
      
      #Update probabilities based on newly updated propensities
      for(j in 1:as.numeric(M[n])){
        p[n,j,t] = q[n,j,t]/sum(q[n,1:as.numeric(M[n]),t])     
      } #next j in the update of probability loop
      
    } #next n
    
  } #next t
  
  #Write the values of p, q, and k to global arrays for later use and analysis
  P[,,,i]=p
  Q[,,,i]=q
  K[,,i]=k
  
} #next i


###Recreate the Plots from Erev and Roth

#Initialize the average arrays
P_avg = array(data = NA, dim = c(NumPlayer,max(as.numeric(M[1:NumPlayer])),Rounds))
Q_avg = P_avg

#Loop through the players and choices and average over the iterations observed at each round
for (n in 1:NumPlayer){
  for (j in 1:as.numeric(M[n])) {
    
    SubP=P[n,j,,]
    P_avg[n,j,]=rowMeans(SubP)
    
    SubQ=Q[n,j,,]
    Q_avg[n,j,]=rowMeans(SubQ)
    
  }
}



###### PLOTTING FUNCTIONS #####

##Plot ALL the iterations
matplot(1:Rounds,P[1,1,,], typ = "l", pch=19, col = "red",
        xlim=c(0,100), ylim=c(0,1),
        xlab = "Round (t)", ylab = "Probability of Choosing Choice 1"
)
matlines(P[2,1,,], pch = 18, col="blue")
legend("topright", legend=c("Player 1", "Player 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       text.font=4)


##Plot a single random iteration
it = sample(1:Iterations, 1)
plot(1:Rounds,P[1,1,,it], typ = "l", pch=19, col = "red",
     xlim=c(0,100), ylim=c(0,1),
     xlab = "Round (t)", ylab = "Probability of Choosing Choice 1"
)
lines(P[2,1,,it], pch = 18, col="blue")
legend("topright", legend=c("Player 1", "Player 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       text.font=4)


##Plot the average of the probabilities for each player

plot(1:Rounds,P_avg[1,1,], typ = "l", pch=19, col = "red", lwd=2,
     xlim=c(0,Rounds), ylim=c(0,1),
     xlab = "Round (t)", ylab = "Probability of Choosing Choice 1"
)
lines(1:Rounds,P_avg[2,1,], pch = 18, col="blue", lty = 2, lwd=2)
legend("topright", legend=c("Player 1", "Player 2"), lwd = 2,
       col=c("red", "blue"), lty=1:2, cex=0.8,
       text.font=4)



##Plot Binned data to match Erev and Roth outputs

#Bin the iterations of data
Bins = 8
Bin_Size = round(Rounds/Bins,0)

t_vector = seq(Rounds/Bins/2, Rounds, Rounds/Bins)
t_vector = append(1, t_vector)

P_bin = array(data = NA, dim = c(NumPlayer,max(as.numeric(M[1:NumPlayer])),Bins+1))
Q_bin = P_bin


#Loop through the players and choices and create an average bin by round
for (n in 1:NumPlayer){
  for (j in 1:as.numeric(M[n])) {
    P_bin[n,j,]=append(P_avg[n,j,1],colMeans(matrix(P_avg[n,j,],Bin_Size)))
    Q_bin[n,j,]=append(Q_avg[n,j,1],colMeans(matrix(P_avg[n,j,],Bin_Size)))
  }
}

#Plot the binned data
plot(t_vector,P_bin[1,1,], typ = "b", pch=19, col = "red", lwd=2,
     xlim=c(0,Rounds), ylim=c(0,1),
     xlab = "Round (t)", ylab = "Probability of Choosing Choice 1"
)
lines(t_vector,P_bin[2,1,], typ = "b", pch = 18, col="blue", lty = 2, lwd=2)
legend("topright", legend=c("Player 1", "Player 2"), lwd = 2,
       col=c("red", "blue"), lty=1:2, cex=0.8,
       text.font=4)



