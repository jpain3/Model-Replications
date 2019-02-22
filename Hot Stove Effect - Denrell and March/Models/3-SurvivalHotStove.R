#Model Replication
#Adaptation as Information Restriction
#Jerker Denrell, James G. March, 2001

#Created by James Paine
#15.879 - Simulation Models in Social and Behavioral Sciences
#February 22, 2019


#####PART 3 - COMPETITIVE SURVIVAL WITHOUT COMPETENCY BUILDING#####

###System Parameters

#Expectation of the 'risky' or 'new' alternative
X = 10
#Standard deviation of the 'risky' or 'new' alternative
S = 10

#Expectation of the 'certain' alternative
Y = 10

###Simulation Paramters
replications = 5
firms = 2
periods = 50
agents = 100
RiskyFract = 0.5

#See page 529 for descriptions of each reproduction mechanism
#  1 = Uniformly random
#  2 = proportional to number of surviving firms
#  3 = proportional to total performance of surviving firms
#  4 = proportional to average performance of surviving firms
ReproductionMechanism = 2

#Set to FALSE to match the paper
#Do the surviving firms keep their previous performance round-by-round?
CummulativePerformance = FALSE

#Display full iteration details - adds time, only use for small replications
fulldetail = 0

#Create arrays to range over for later plotting
W = seq(from = .05, to = 0.95, by = ((0.95-0.05)/18))
H = c(0.3,0.5,1)


#t = time period (note, the paper is indexed relative to t=0, but here it's relative to t=1)
#w = fraction of the population eliminated each time period
#h = affects the sensitivity of reproduction to the aggregate performance of a type (see page 529)

#Create a data frame to keep track of the performance of each firm
PerfList = data.frame(matrix(NA, ncol =4, nrow = agents), stringsAsFactors = FALSE)
colnames(PerfList) = c("Agent","Type","Performance","Rank")

#Create an array to keep track of the fraction of X (or 1) type firms at the end of each run
FractX = array(data=NA, dim=c(length(W),length(H), replications))

t=1
n = 1
w = 0.3
h = 0.3
sens = 1
elim = 1

starttime = proc.time()

for (rep in 1:replications) {
  
  if (fulldetail != 1) {
    cat("\014")
    print(paste0("Replication ",rep," of ", replications))
    print(paste0("Elapsed time since last replication: ",(proc.time() - starttime)[3]))
    print(paste0("Avg Replication time: ",(proc.time() - starttime)[3]/rep," per replication"))
    print(paste0("Est Time Remaining: ",((proc.time() - starttime)[3]/rep)*(replications-rep)))
    
    
  }
  
  for (sens in 1: length(H)) {
  
    h = H[sens]
    
    for (elim in 1:length(W)) { #step through each elimination percentage w
    
      w = W[elim]
      
      #Initialize the data matrix for firm performance of each firm type
      PerfList$Agent = seq.int(1,agents)
      PerfList$Performance = 0
      PerfList$Type = append(rep(1,round(RiskyFract*agents)),rep(2,agents-round(RiskyFract*agents)))
      
      #Get the number of frims to elminate each round based on the value of w
      kills = round(w*nrow(PerfList))
      
      for (t in 1:(periods+1)) {
        
        if (fulldetail == 1) {
          cat("\014")
          print(paste0("Replication ",rep," of ", replications))
          print(paste0("h value: ", h))
          print(paste0("w value: ", w))
          print(paste0("Time Period: ", (t-1), " of ", periods))
        }  
          
        for (n in 1:agents) {
          
          #Determine the firm performance based on its type
          
          
          if (CummulativePerformance == TRUE) {
            CummPerf = 1  
          } else{
            CummPerf = 0
          }
          
          if (PerfList$Type[n] == 1) {
            PerfList$Performance[n] = (CummPerf*PerfList$Performance[n]) + rnorm(1, mean = X, sd = S)
          }
          
          if (PerfList$Type[n] == 2) {
            PerfList$Performance[n] = (CummPerf*PerfList$Performance[n]) + Y
          }
          
        } #next agent n
       
        #Order based on performance from worst to best
        PerfList$Rank = rank(PerfList$"Performance",ties.method = "random")
        PerfList = PerfList[order(PerfList$"Rank"),]
        
        SurviveList = PerfList[(kills+1):agents,]
      
        #Determine the various performance factors for the survivors that affect reproduction
        
        T1 = sum((SurviveList$Type == 1)*SurviveList$Performance)
        T2 = sum((SurviveList$Type == 2)*SurviveList$Performance)
        
        N1 = sum((SurviveList$Type == 1))
        N2 = sum((SurviveList$Type == 2))
        
        t
        
        A1 = T1/N1
        A2 = T2/N2
        
        #Define reproduction probability based on user choices
        
        if (ReproductionMechanism == 1) { #uniformly random replacement
          r1 = 1/firms
        } else {  #replacements proportional to the ammount/performance of firms
          
          #Avoid erroneous rates when one population is totally eliminated
          if (N1 == 0) {
            r1 = 0
          } else if (N2 == 0) {
            r1 = 1
          } else {
            
            if (ReproductionMechanism == 2) {
              r1 = N1^h/(N1^h+N2^h)
            } else if (ReproductionMechanism ==3) {
              r1 = T1^h/(T1^h+T2^h)   
            } else if (ReproductionMechanism ==4) {
              r1 = A1^h/(A1^h+A2^h)
            }
            
          }
        } 
        
        #Generate new firms based on the reproduction probabilities  
        NewTypes = sample(c(1:2),kills,replace=TRUE,prob=c(r1,1-r1))
        
        #Record the new firm types
        PerfList[1:kills,]$Type = NewTypes
        #Reset the new firms performance to 0
        PerfList[1:kills,]$Performance = 0
        
      } # Next time period t
      
      FractX[elim,sens,rep] = sum(PerfList$Type == 1)/agents
      
    
    } # next elimination percentage w  
  
  } # next sensitivity factor h

} # next replication

endtime = proc.time()

ElaspedTime = endtime-starttime
ElaspedTime

#Average each replication and determine the standard deviation
FractX_avg = apply(FractX, c(1,2), mean)
FractX_sd = apply(FractX, c(1,2), sd)

rep=1
matplot(W,FractX[,,rep], typ = "l", pch=4, col = "black", lwd=1, lty =1,
        xlim=c(0,1), ylim=c(0,1),
        xlab = "Percent Eliminated Each Round (w)", ylab = "Fraction of Risk-Choosing Agents at Round 50"
)


matplot(W,FractX_avg, typ = "l", pch=4, col = "black", lwd=1, lty =1,
        xlim=c(0,1), ylim=c(0,1),
        xlab = "Percent Eliminated Each Round (w)", ylab = "Fraction of Risk-Choosing Agents at Round 50"
)


##Plot of Figure 3 on page 530

plot(W,FractX_avg[,1], typ = "l", col = "black", lwd=1, lty =1,
        xlim=c(0,1), ylim=c(0,1),
        xlab = "Percent Eliminated Each Round (w)", ylab = "Fraction of Risk-Choosing Agents at Round 50")
    points(W,FractX_avg[,1], typ = "p", pch=23, lwd=1, bg = "white")
    lines(W,FractX_avg[,2], type = "l", col = "black", lwd=1, lty = 3)
    points(W,FractX_avg[,2], typ = "p", pch=8, lwd=1.5)
    lines(W,FractX_avg[,3], type = "l", col = "black", lwd=1, lty = 3)
    points(W,FractX_avg[,3], typ = "p", pch=21, lwd=1.5, bg = "white")
    #text(c(0.25,0.38,0.5),c(0.28,0.13,0.08),paste("h = ", H), cex = .75)
    
    
