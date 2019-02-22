#Model Replication
#Adaptation as Information Restriction
#Jerker Denrell, James G. March, 2001

#Created by James Paine
#15.879 - Simulation Models in Social and Behavioral Sciences
#February 22, 2019

#####PART 4 - COMPETITIVE SURVIVAL WITH COMPETENCY BUILDING#####

###System Parameters

#Expectation of the 'risky' or 'new' alternative
X = 15
#Standard deviation of the 'risky' or 'new' alternative
S = 5

#Expectation of the 'certain' alternative
Y = 10

###Simulation Paramters
replications = 5
firms = 2
periods = 50
agents = 100
RiskyFract = 0.5
h = 0.5
k = 0.5
c0 = 0.3

#See page 529 for descriptions of each reproduction mechanism
#  1 = Uniformly random
#  2 = proportional to number of surviving firms
#  3 = proportional to total performance of surviving firms
#  4 = proportional to average performance of surviving firms
ReproductionMechanism = 2

#Set to FALSE to match the paper
#Do the surviving firms keep their previous performance round-by-round?
CummulativePerformance = TRUE

LearningTransfer = FALSE

#Display full iteration details - adds time, only use for small replications
fulldetail = 0

#Create arrays to range over for later plotting
W = seq(from = .05, to = 0.95, by = ((0.95-0.05)/18))
#D = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
D = c(0.1,0.4,0.6)

#t = time period (note, the paper is indexed relative to t=0, but here it's relative to t=1)
#w = fraction of the population eliminated each time period
#h = affects the sensitivity of reproduction to the aggregate performance of a type (see page 529)

#Create a data frame to keep track of the performance of each firm
PerfList = data.frame(matrix(NA, ncol =3, nrow = agents), stringsAsFactors = FALSE)
colnames(PerfList) = c("Type","Performance","c")

#Create an array to keep track of the fraction of X (or 1) type firms at the end of each run
FractX = array(data=NA, dim=c(length(W),length(D), replications))

starttime = proc.time()

for (rplct in 1:replications) {
  
  if (fulldetail != 1) {
    cat("\014")
    print(paste0("Replication ",rplct," of ", replications))
    print(paste0("Elapsed time since last replication: ",(proc.time() - starttime)[3],"s"))
    print(paste0("Avg Replication time: ",(proc.time() - starttime)[3]/rplct,"s per replication"))
    print(paste0("Est Time Remaining: ",((proc.time() - starttime)[3]/rplct)*(replications-rplct),"s"))
    
    
  }
  
  
  for (LearnParm in 1:length(D)) {
  
    d = D[LearnParm]

    for (elim in 1:length(W)) { #step through each elimination percentage w

      w = W[elim]
     
      #Initialize the data matrix for firm performance of each firm type
      PerfList$Performance = 0
      PerfList$Type = append(rep(1,round(RiskyFract*agents)),rep(2,agents-round(RiskyFract*agents)))
      
      #Set the initial competency to c0 before stepping forward through time periods
      PerfList$c = append(rep(c0,round(RiskyFract*agents)),rep(NA,agents-round(RiskyFract*agents)))
      
      #Get the number of frims to elminate each round based on the value of w
      kills = round(w*nrow(PerfList))
      
      for (t in 1:(periods)) {

        if (fulldetail == 1) {
          cat("\014")
          print(paste0("Replication ",rplct," of ", replications))
          print(paste0("d value: ", d))
          print(paste0("w value: ", w))
          print(paste0("Time Period: ", (t-1), " of ", periods))
        }  

        for (n in 1:agents) {
          
          if (CummulativePerformance == TRUE) {
            CummPerf = 1  
          } else{
            CummPerf = 0
          }
          
          
          #Determine the firm performance based on its type
          
          if (PerfList$Type[n] == 1) {
            #Get the st dev based on the current competency c
            stdev = (S/PerfList$c[n])^k
            avg = PerfList$c[n]*X
            #store the performance for this agent
            PerfList$Performance[n] = (CummPerf*PerfList$Performance[n]) + rnorm(1, mean = avg, sd = stdev)
            
            #Update the agent's competency with the risky process for the use in the next round (if they survive)
            PerfList$c[n] = PerfList$c[n] + d*(1-PerfList$c[n])
            
          }
          
          if (PerfList$Type[n] == 2) {
            
            #For these agents, the performance is a constant
            PerfList$Performance[n] = (CummPerf*PerfList$Performance[n]) + Y

          }
          
        } #next agent n
       
        #Order, from worst to best, based on performance from worst to best
        #PerfList$Rank = rank(PerfList$"Performance",ties.method = "random")
        #PerfList = PerfList[order(PerfList$"Rank"),]
        PerfList = PerfList[order(PerfList$"Performance"),]
        
        #Get list of suriving firms (bottom of ordered list)
        SurviveList = PerfList[(kills+1):agents,]
        
        #Determine the various performance factors for the survivors that affect reproduction
        
        T1 = sum((SurviveList$Type == 1)*SurviveList$Performance)
        T2 = sum((SurviveList$Type == 2)*SurviveList$Performance)
        
        N1 = sum((SurviveList$Type == 1))
        N2 = sum((SurviveList$Type == 2))

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
        #Reset the new firms performance to NA
        PerfList[1:kills,]$Performance = 0

        #Reset new firms of type 1 to the baseline compentency with the process
        
        AvgC = mean(SurviveList[SurviveList$Type==1,]$c)
        
        if (LearningTransfer == TRUE) {
          NewC = AvgC
        } else {
          NewC = c0
        }
        
        NewComp = replace(replace(NewTypes, NewTypes==2, NA),NewTypes==1, NewC)
        PerfList[1:kills,]$c = NewComp
        
        sum(PerfList$Type == 1)/agents
        
      } # Next time period t
      
      FractX[elim,LearnParm,rplct] = sum(PerfList$Type == 1)/agents


    } # next elimination percentage w  
  
  } # next learning speed paramter d

} # next replication

endtime = proc.time()

(ElaspedTime = endtime-starttime)

#Average each replication and determine the standard deviation
FractX_avg = apply(FractX, c(1,2), mean)
FractX_sd = apply(FractX, c(1,2), sd)

##Plot of Figure 4 on page 530


matplot(W,FractX_avg[,], typ = "l", col = "black", lwd=1, lty =1,
     xlim=c(0,1), ylim=c(0,1),
     xlab = "Percent Eliminated Each Round (w)", ylab = "Fraction of Risk-Choosing Agents at Round 50")

##Figure 4

plot(W,FractX_avg[,1], typ = "l", col = "black", lwd=1, lty =1,
        xlim=c(0,1), ylim=c(0,1),
        xlab = "Percent Eliminated Each Round (w)", ylab = "Fraction of Risk-Choosing Agents at Round 50")
    lines(W,FractX_avg[,2], type = "l", col = "black", lwd=1, lty = 1)
    lines(W,FractX_avg[,3], type = "l", col = "black", lwd=1, lty = 3)
    lines(W,FractX_avg[,4], type = "l", col = "black", lwd=1, lty = 4)
    lines(W,FractX_avg[,5], type = "l", col = "black", lwd=1, lty = 6)
    lines(W,FractX_avg[,6], type = "l", col = "black", lwd=1, lty = 1)
    lines(W,FractX_avg[,7], type = "l", col = "black", lwd=1, lty = 6)
    lines(W,FractX_avg[,8], type = "l", col = "black", lwd=1, lty = 6)
    lines(W,FractX_avg[,9], type = "l", col = "black", lwd=1, lty = 6)
    lines(W,FractX_avg[,10], type = "l", col = "black", lwd=1, lty = 1)
    points(W,FractX_avg[,1], typ = "p", pch=23, lwd=1, bg = "white")
    points(W,FractX_avg[,2], typ = "p", pch=21, lwd=1.5, bg = "black")
    points(W,FractX_avg[,3], typ = "p", pch=8, lwd=1.5)
    points(W,FractX_avg[,4], typ = "p", pch=25, lwd=1.5, bg = "white")
    points(W,FractX_avg[,5], typ = "p", pch=10, lwd=1.5)
    points(W,FractX_avg[,6], typ = "p", pch=4, lwd=1.5)
    text(c(0.9,0.8,0.65,.53,.4,.34,.22),c(0.05,0.1,0.15,.3,.4,.55,.8),paste("d = ", D[c(1:6,10)]), cex = .75)
    
