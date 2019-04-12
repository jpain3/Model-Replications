#James Paine
#NK Model
#Adaptation on Rugged Landscapes
#Levinthal 1997

#Load library for cominatorics
library(gtools)

set.seed(123)

#Number of attributes
N = 10

Horizon = 50

#Number of entities to simulate
NumEntity = 100

#Note that each attribute can have a value of 0 or 1

#Degree of interaction between attributues and their neighbors
K_values=c(0,1,5)

replications = 100

start_time <- Sys.time()

#Initialize the vector of unique firms
Unique_Firms = matrix(100,nrow = (Horizon+1),ncol = replications)

Unique_Firms_avg = matrix(100,nrow = (Horizon+1),ncol = length(K_values))
colnames(Unique_Firms_avg) = K_values

q = 1

for (q in 1:length(K_values)) {
  
  K = K_values[q] 
  
  #Number of possible fitness values that could occur for each element given the subsequent element values
  NumFitness = 2^(K+1)
  
  print(paste("K value =",K))
  
  FullTable = data.frame(permutations(n=2,r=N,v=c(0,1), repeats.allowed=T))
  FullTable$MatchString = apply(FullTable[,1:N],1,paste, collapse=" ")
  
  r =1
  
  for (r in 1:replications) {
    
    
    print(paste("K value =",K,"Replication =",(r)))
    
    print("Determining Permuation Values...")
    
    #Generate a permuation table with corresponding fitness values based on K
    FitTable = data.frame(permutations(n=2,r=(K+1),v=c(0,1), repeats.allowed=T))
    FitTable$MatchString = apply(FitTable,1,paste, collapse=" ")
    FitTable$Fitness = runif(NumFitness)
    
    #Determine the fitness for every possible entity
    
    print("Populating Fitness Lookup Table...")
    FullTable$Fitness = rep(0,nrow(FullTable))
    
    for (n in 1:nrow(FullTable)) {
      
      EntityString = unname(unlist(FullTable[n,]))
      
      #Calculate the total fitness of the entity
      TotalFitness = 0
      
      #Loop through each attribute
      for (i in 1:N) {
        
        TestString = ""
        
        #Loop over the number of interaction points
        for (k in 0:K) {
          
          #Allow for looping when at the end of the string
          index = max(i+k - N, 0)
          
          if (index>0) {
            index = index
          } else {
            index = i+k
          }
          
          TestString = trimws(paste(TestString,EntityString[index],sep=" "))
        } 
        
        #Lookup the fitness value for this string
        CurrentFitness = FitTable[(FitTable$MatchString == TestString),]$Fitness
        TotalFitness = CurrentFitness+TotalFitness
        
      }
      
      #Write the total Fitness to the Entity Matrix
      FullTable[n,"Fitness"] = TotalFitness/N
      
    }  
    
    #Look up nearest neighbors and the one with the maximum fitness value
    print(paste("Generating list of maximum fitness neighbors..."))
    Max_Neighbor = function(EntityString) {
      DiffMatrix = t(abs(apply(FullTable[,1:N],1,'-',unname(unlist(EntityString)))))
      DiffVector = apply(DiffMatrix,1,sum)
      SubsetDiffFitness = which(DiffVector == 1)
      NeighborFitness = FullTable[SubsetDiffFitness,"Fitness"]
      MaxNeighborIndex = SubsetDiffFitness[which.max(NeighborFitness)]
      return(MaxNeighborIndex)
    }
    
    MaxNeighborLookup = apply(FullTable[,1:N],1,Max_Neighbor)
    
    
    #Initialize the matrix of entities with random attributes
    EntityMatrix = data.frame(matrix(sample(c(0,1), replace=TRUE, size = N*NumEntity),NumEntity))
    EntityMatrix$MatchString = apply(EntityMatrix[,1:N],1,paste, collapse=" ")
    
    #Match the fitness values to the entities from the full look-up table
    EntityMatrix$Fitness = FullTable[match(EntityMatrix$MatchString,FullTable$MatchString),"Fitness"]
    
    #Initial number Number of unique entries in the table:
    Unique_Firms[1,r] = dim(table(EntityMatrix$MatchString))
    #Unique_Firms[1,r] = dim(table(EntityMatrix$Fitness))
    
    print(paste("Iterating over ",Horizon,"timesteps..."))
    
    for (t in 2:(Horizon+1)) {
      
      #print(paste("time step =",(t-1)))
      
      #Determine total population Genetic Load:
      GL = 1 - mean(EntityMatrix$Fitness)/max(EntityMatrix$Fitness)
      
      #Determine the bottom performers to cull based on Genetic Load
      DeathToll = sum(sample(c(0,1), replace=TRUE, size = NumEntity, prob = c(1-GL,GL)))
      
      if (DeathToll > 0) {
        NumSurvive = NumEntity - DeathToll
        DeathList = order(EntityMatrix$Fitness)[1:DeathToll]
  
        #Cull the Entity Matrix to just the survivers
        EntityMatrix = EntityMatrix[-DeathList,]
        rownames(EntityMatrix) = seq(from=1,to=nrow(EntityMatrix),by=1)
      }
        
      ########
      
      #LOCAL AND LONG DISTANCE SEARCH
      
      #######

      
      if (DeathToll >0) {
        
        ######
        #REFILL ENTITY TABLE
        ######
        
        #Determine surviving population Genetic Load:
        GL = 1 - mean(EntityMatrix$Fitness)/max(EntityMatrix$Fitness)
        
        #Determine propability of replication of survivors
        Prob_Rep = EntityMatrix$Fitness/sum(EntityMatrix$Fitness)
        
        #Get number of new organizatons that are birthed via replication (versus randomization):
        Num_Replicants = sum(sample(c(1,0), replace=TRUE, size = DeathToll, prob = c(1-GL,GL)))
        Num_Randos = DeathToll - Num_Replicants
        
        #Get the indicies of the entities to replicate and create a matrix of these replicated entries
        if (Num_Replicants > 0) {
          Replicants = sample(seq(from=1,to=NumSurvive,by=1), replace=TRUE, size = Num_Replicants, prob = Prob_Rep)
          ReplicantMatrix = EntityMatrix[Replicants,]
        } else {
          ReplicantMatrix = EntityMatrix[FALSE,]
        }
  
        ##Create the needed number of random new entities
        if (Num_Randos > 0) {
          #Initialize the matrix of entities with random attributes
          RandoMatrix = data.frame(matrix(sample(c(0,1), replace=TRUE, size = N*Num_Randos),Num_Randos))
          RandoMatrix$MatchString = apply(RandoMatrix[,1:N],1,paste, collapse=" ")
          #Match the fitness values to the entities from the full look-up table
          RandoMatrix$Fitness = FullTable[match(RandoMatrix$MatchString,FullTable$MatchString),"Fitness"]
        } else {
          RandoMatrix = EntityMatrix[FALSE,]
        }
          
      
        #Combine Surivors, Replicants, and Randos back together
        EntityMatrix = rbind(EntityMatrix,ReplicantMatrix,RandoMatrix)
        rownames(EntityMatrix) = seq(from=1,to=nrow(EntityMatrix),by=1)
      
      }
      
      #Number of unique entries in the table:
      Unique_Firms[t,r] = dim(table(EntityMatrix$MatchString))
      
    } #next t
    
    print(paste(""))
    
  } #next r
  
  Unique_Firms_avg[,q] = apply(Unique_Firms, 1, mean)
  
} #next k_value


end_time <- Sys.time()

print(end_time - start_time)



plot((0:Horizon),Unique_Firms_avg[,3], typ = "l", pch=4, col = "black", lwd=2, lty =1,
     xlim=c(0,Horizon), ylim=c(0,NumEntity),
     xlab = "Time", ylab = "Organizational Forms"
)


matplot((0:Horizon),Unique_Firms_avg, typ = "l", pch=4, col = "grey", lwd=1, lty =1,
        xlim=c(0,Horizon), ylim=c(0,NumEntity),
        xlab = "Time", ylab = "Organizational Forms"
)

plot((0:Horizon),Unique_Firms_avg[,1], typ = "l", pch=19, col = "black", lwd=1,
     xlim=c(0,Horizon), ylim=c(0,NumEntity),
     xlab = "Time", ylab = "",
     main = "Organizational Forms"
)
lines((0:Horizon),Unique_Firms_avg[,2], typ = "l", pch = 18, col="black", lty = 2, lwd=1)
lines((0:Horizon),Unique_Firms_avg[,3], typ = "l", pch = 18, col="black", lty = 3, lwd=1)
legend("topright", legend=paste("k =",K_values), lwd = 1,
       col="black", lty=1:3, cex=0.8,
       text.font=4)
