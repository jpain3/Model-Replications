#James Paine
#NK Model
#Adaptation on Rugged Landscapes
#Levinthal 1997

#Load library for cominatorics
library(gtools)

#set.seed(123)

#Number of attributes
N = 10

Horizon = 10

#Number of entities to simulate
NumEntity = 100

#Note that each attribute can have a value of 0 or 1

#Degree of interaction between attributues and their neighbors
K_values=c(1)

replications = 1

start_time <- Sys.time()

#Initialize the vector of unique firms
Unique_Firms = matrix(100,nrow = (Horizon+1),ncol = replications)

Unique_Firms_avg = matrix(100,nrow = (Horizon+1),ncol = length(K_values))
colnames(Unique_Firms_avg) = K_values

repeat{

  for (q in 1:length(K_values)) {
  
    K = K_values[q] 
    
    #Number of possible fitness values that could occur for each element given the subsequent element values
    NumFitness = 2^(K+1)
    
    print(paste("K value =",K))
    print("Determining Permuation Values...")
    
    #Generate a permuation table with corresponding fitness values based on K
    FitTable = data.frame(permutations(n=2,r=(K+1),v=c(0,1), repeats.allowed=T))
    FitTable$MatchString = apply(FitTable,1,paste, collapse=" ")
    FitTable$Fitness = runif(NumFitness)
    
    
    FullTable = data.frame(permutations(n=2,r=N,v=c(0,1), repeats.allowed=T))
    FullTable$MatchString = apply(FullTable[,1:N],1,paste, collapse=" ")
    FullTable$Fitness = rep(0,nrow(FullTable))
    
    #Determine the fitness for every possible entity
    
    print("Populating Fitness Lookup Table...")
    
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
    
    
    for (r in 1:replications) {
    
      
      print(paste("Replication =",(r)))
      
      #Initialize the matrix of entities with random attributes
      EntityMatrix = data.frame(matrix(sample(c(0,1), replace=TRUE, size = N*NumEntity),NumEntity))
      EntityMatrix$MatchString = apply(EntityMatrix[,1:N],1,paste, collapse=" ")
      
      #Match the fitness values to the entities from the full look-up table
      EntityMatrix$Fitness = FullTable[match(EntityMatrix$MatchString,FullTable$MatchString),"Fitness"]
    
      #Initial number Number of unique entries in the table:
      Unique_Firms[1,r] = dim(table(EntityMatrix$MatchString))
      #Unique_Firms[1,r] = dim(table(EntityMatrix$Fitness))
      
      for (t in 2:(Horizon+1)) {
    
        print(paste("time step =",(t-1)))
    
        #Local Scanning
        
        n=1
    
        for (n in 1:NumEntity) {
          
          #print(paste("Entity =",n))
          
          #Get the entity that is looking for local neighbors
          EntityString = unname(unlist(EntityMatrix[n,1:N]))
        
          #Determine how different all other entities are from current entity
          DiffMatrix = t(abs(apply(FullTable[,1:N],1,'-',EntityString)))
          DiffVector = apply(DiffMatrix,1,sum)
          
          #Create a data frame containing the degree of difference and the fitness of all entities
          DiffFitnessMatrix = data.frame(cbind(DiffVector,FullTable[,"Fitness"]))
          colnames(DiffFitnessMatrix) = c("DiffVector","Fitness")
          
          #Find the subset that are only one degree off
          SubsetDiffFitness = DiffFitnessMatrix[DiffFitnessMatrix$DiffVector==1,]
          
          #Find the index of first entity that has greater fitness
          FirstLargest = match(TRUE,(SubsetDiffFitness$Fitness > EntityMatrix[n,"Fitness"]))
          
          #If NA, keep the entity the same, else swap it out
          if (is.na(FirstLargest) == TRUE) {
            EntityMatrix[n,] = EntityMatrix[n,]
          } else {
            SwapRow = as.integer(rownames(SubsetDiffFitness[FirstLargest,]))
            #Make the swap
            EntityMatrix[n,] = FullTable[SwapRow,]
            
          }
          
    
          
        }
      
        #Number of unique entries in the table:
        Unique_Firms[t,r] = dim(table(EntityMatrix$MatchString))
        #Unique_Firms[t,r] = dim(table(EntityMatrix$Fitness))
        #write.table(apply(EntityMatrix[,1:10],1,paste, collapse=" "), "clipboard", sep="\t") #Export to Clipboard for copy/paste into external docs
        
      } #next t
      
    } #next r
    
    Unique_Firms_avg[,q] = apply(Unique_Firms, 1, mean)
    
    
    
  } #next k_value

  print(paste("Unique Fitnesses in Existance =",dim(table(EntityMatrix$Fitness))))
  
  if (dim(table(EntityMatrix$Fitness))>2) {
    
    break
    
  }
  
} #loop
  
end_time <- Sys.time()

print(end_time - start_time)

BarChartDataFrame = data.frame(table(EntityMatrix$Fitness))
BarChartData = BarChartDataFrame[,2]
names(BarChartData)=round(as.numeric(as.vector(BarChartDataFrame[,1])),3)

barplot(BarChartData, main="Number of Organizations", 
        xlab="Fitness", ylim = c(0,100))

