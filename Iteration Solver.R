source("General attack solver.R")

IterationSolver<-function(Omega,AttackTimeDistribution,B,b,Cost,Lambda,TypeOfAttackTimeDis="CDF")
{
 if(TypeOfAttackTimeDis=="PDF")
 {
   #Conversion is required
   AttackTimeDistribution=GenerateCDF(AttackTimePDFDistribution)
 }
  
 #Generate cost matrix
 CostToProgressMatrix=GenerateCostToProgressMatrix(B+1,b+1,Cost,Lambda,AttackTimeDistribution)
 
 #Then we peform iterations to find g with an initial guess of c * lambda
 InitialGuess= Cost * Lambda
 EvaluatedEquilibrium=InitialGuess
 MinTolerance=0.0001
 MaxSteps=1000
 Step=1
 Tolerance=MinTolerance+1
 CurrentPlan=FindRenewInMatrix(InitialGuess,CostToProgressMatrix)
 while(Tolerance>MinTolerance && Step<MaxSteps)
 {
   OldEquilibriumValue=EvaluatedEquilibrium
   EvaluatedEquilibrium=FindEquilibriumValue(CurrentPlan,CostToProgressMatrix,Lambda,Omega)
   CurrentPlan=FindRenewInMatrix(EvaluatedEquilibrium,CostToProgressMatrix)
   
   #Now see if equilibrium value has changed (Equally we could see if plan has changed ??)
   
   Tolerance=abs(EvaluatedEquilibrium-OldEquilibriumValue)/OldEquilibriumValue
   Step=Step+1
 }
 return(list(Plan=CurrentPlan,EquilibriumValue=EvaluatedEquilibrium))
}