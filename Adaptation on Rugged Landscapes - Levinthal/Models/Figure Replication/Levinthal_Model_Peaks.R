#James Paine
#NK Model
#Adaptation on Rugged Landscapes
#Levinthal 1997

#Load library for cominatorics
library(gtools)

set.seed(123)

#Number of attributes
N = 10

#Note that each attribute can have a value of 0 or 1

#Degree of interaction between attributues and their neighbors
K_values=seq(from = 0, to = (N-1), by = 1)


replications = 100

start_time <- Sys.time()

Peak_Counts = matrix(0,nrow = length(K_values),ncol = replications)


for (q in 1:length(K_values)) {

  K = K_values[q] 
  
  #Number of possible fitness values that could occur for each element given the subsequent element values
  NumFitness = 2^(K+1)
  
  print(paste("K value =",K))
  
  FullTable = data.frame(permutations(n=2,r=N,v=c(0,1), repeats.allowed=T))
  FullTable$MatchString = apply(FullTable[,1:N],1,paste, collapse=" ")
  
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

    #Get count of local peaks
    Local_Peaks = sum((FullTable$Fitness > FullTable[MaxNeighborLookup,"Fitness"]), na.rm = TRUE)
    
    print(paste("Local Peaks Found = ",Local_Peaks))
    print(paste(""))
    
    Peak_Counts[q,r] = Local_Peaks

    
  } #next r
  

} #next k_value


end_time <- Sys.time()

print(end_time - start_time)

Peak_Counts_Avg = data.frame(Avg = rep(NA,length(K_values)))
Peak_Counts_Avg$Avg = apply(Peak_Counts, 1, mean)
Peak_Counts_Avg$SD = apply(Peak_Counts, 1, sd)
Peak_Counts_Avg$"Avg+SD" = Peak_Counts_Avg$Avg + Peak_Counts_Avg$SD
Peak_Counts_Avg$"Avg-SD" = Peak_Counts_Avg$Avg - Peak_Counts_Avg$SD

plot(K_values,Peak_Counts_Avg$Avg, typ = "l", pch=19, col = "black", lwd=1,
     xlim=c(0,1.1*max(K_values)), ylim=c(0,1.1*max(Peak_Counts_Avg$Avg,Peak_Counts_Avg$"Avg+SD")),
     xlab = "K Value", ylab = "",
     main = "Number of Peaks"
)
lines(K_values,Peak_Counts_Avg$"Avg+SD", typ = "l", pch = 18, col="black", lty = 2, lwd=1)
lines(K_values,Peak_Counts_Avg$"Avg-SD", typ = "l", pch = 18, col="black", lty = 3, lwd=1)
legend("bottomright", legend=c("Mean","Mean+St.Dev","Mean-St.Dev"), lwd = 1,
       col="black", lty=1:3, cex=0.8,
       text.font=4)

