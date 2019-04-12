#James Paine
#NK Model
#Adaptation on Rugged Landscapes
#Levinthal 1997

#Load library for cominatorics
library(gtools)

#set.seed(10261986)
set.seed(123)

#Number of attributes
N = 10

Horizon = 50

#Number of entities to simulate
NumEntity = 100

#Note that each attribute can have a value of 0 or 1

#Degree of interaction between attributues and their neighbors
K_values=c(1)

OrgAdapation = 0

replications = 1

start_time <- Sys.time()

#Initialize the vector of unique firms
Unique_Firms = matrix(100,nrow = (Horizon+1),ncol = replications)

Unique_Firms_avg = matrix(100,nrow = (Horizon+1),ncol = length(K_values))
colnames(Unique_Firms_avg) = K_values


#Initialize the table of fitness counts

for (q in 1:length(K_values)) {
  
  K = 1
  
  K = K_values[q] 
  
  #Number of possible fitness values that could occur for each element given the subsequent element values
  NumFitness = 2^(K+1)
  
  print(paste("K value =",K))
  
  FullTable = data.frame(permutations(n=2,r=N,v=c(0,1), repeats.allowed=T))
  FullTable$MatchString = apply(FullTable[,1:N],1,paste, collapse=" ")
  
  #Note, in the below the first index is the Entity index in the FullTable, and the second index is time
  DistributionArray = array(dim = c(2^N,Horizon+1,replications))
  
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
    
    #Get the count of Entities in the FullTable for this time step:
    IndexOfEntities = as.integer(rownames(table(match(EntityMatrix$MatchString,FullTable$MatchString))))
    CountOfEntities = unname(table(match(EntityMatrix$MatchString,FullTable$MatchString)))
    DistributionArray[IndexOfEntities,1,r] = CountOfEntities
    
    #Initial number Number of unique entries in the table:
    Unique_Firms[1,r] = dim(table(EntityMatrix$MatchString))
    #Unique_Firms[1,r] = dim(table(EntityMatrix$Fitness))
    

    print(paste("Iterating over ",Horizon,"timesteps..."))
    
    for (t in 2:(Horizon+1)) {
      
      #print(paste("    Culling the weak..."))
      
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
      
      if (OrgAdapation == 1) {
        
        #LOCAL AND LONG DISTANCE SEARCH
        
        #print(paste("    Local and long distance searching..."))
        
        for (n in 1:nrow(EntityMatrix)) {
          
          #print(paste("Entity =",n))
          
          
          #Get the entity that is looking for local neighbors
          EntityString = unname(unlist(EntityMatrix[n,1:N]))
          
          ##Long Distance Scanning
          
          #Generage Random string to compare to:
          RandEntity = data.frame(t(as.matrix(sample(c(0,1), replace=TRUE, size = N))))
          RandEntity$MatchString = apply(RandEntity[,1:N],1,paste, collapse=" ")
          
          #Look up the fitness for this new Random Entity
          RandEntity$Fitness = FullTable[match(RandEntity$MatchString,FullTable$MatchString),"Fitness"]
          
          #If superior, then update the Entity
          if (RandEntity$Fitness > EntityMatrix[n,"Fitness"] ) {
            
            EntityMatrix[n,] = RandEntity
            
          } else {  #Perform local scanning if random jump failed
            
            ##Local Scanning
            #Get largest fitness of neighbors from lookup
            FullTableIndex = which(FullTable$MatchString == EntityMatrix[n,"MatchString"])
            FullTableNeighbor = MaxNeighborLookup[FullTableIndex]
            
            if (FullTable[FullTableNeighbor,"Fitness"] > FullTable[FullTableIndex,"Fitness"]) {
              #Make the swap
              EntityMatrix[n,] = FullTable[FullTableNeighbor,]
            }
            
          }
          
        }
        
      }
      
      
      
      
      #print(paste("    Reproducing based on fitness"))
      
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
    
      #Get the count of Entities in the FullTable for this time step:
      IndexOfEntities = as.integer(rownames(table(match(EntityMatrix$MatchString,FullTable$MatchString))))
      CountOfEntities = unname(table(match(EntityMatrix$MatchString,FullTable$MatchString)))
      DistributionArray[IndexOfEntities,t,r] = CountOfEntities    
      
    } #next t
    
    
      
    print(paste(""))
    
  } #next r
  
  Unique_Firms_avg[,q] = apply(Unique_Firms, 1, mean)
  
} #next k_value


end_time <- Sys.time()
EllaspedRunTime = end_time - start_time

print(EllaspedRunTime)




#Average over the replications
DistributionArray_Avg = apply(DistributionArray, c(1,2), mean, na.rm = TRUE)
DistributionArray_Avg[is.nan(DistributionArray_Avg)] = NA
max(DistributionArray_Avg[,10], na.rm = TRUE)

#
T_index = 10

InstanceCount = DistributionArray_Avg[,T_index][!is.na(DistributionArray_Avg[,T_index])]
InstanceFitness = FullTable$Fitness[!is.na(DistributionArray_Avg[,T_index])]
TotalFitnessList = rep(InstanceFitness, times = InstanceCount)
df10 = data.frame(table(TotalFitnessList))
df10$t = T_index

T_index = 30

InstanceCount = DistributionArray_Avg[,T_index][!is.na(DistributionArray_Avg[,T_index])]
InstanceFitness = FullTable$Fitness[!is.na(DistributionArray_Avg[,T_index])]
TotalFitnessList = rep(InstanceFitness, times = InstanceCount)
df30 = data.frame(table(TotalFitnessList))
df30$t = T_index

T_index = 50

InstanceCount = DistributionArray_Avg[,T_index][!is.na(DistributionArray_Avg[,T_index])]
InstanceFitness = FullTable$Fitness[!is.na(DistributionArray_Avg[,T_index])]
TotalFitnessList = rep(InstanceFitness, times = InstanceCount)
df50 = data.frame(table(TotalFitnessList))
df50$t = T_index

FullDistDF = merge(merge(df10,df30,all=TRUE),df50,all=TRUE)
colnames(FullDistDF) = c("Fitness","Count","t")



counts = reshape(FullDistDF, idvar = "Fitness", timevar = "t", direction = "wide")
colnames(counts) = gsub("Count.", "T = ", colnames(counts))
counts = counts[ , order(names(counts))]

LegendLabels = round(as.numeric(as.character(counts[,1])),3)

datacounts = as.matrix(counts[,2:4])

barplot(datacounts, beside=TRUE,
        main = "Distribution of Forms (K=1)",
        xlab = "Time",
        ylab = "Number",
        legend = LegendLabels, args.legend = list(x = "top", ncol = 3, bty = "n")
        )

