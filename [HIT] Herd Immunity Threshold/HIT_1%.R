# This code aims to predict when HIT would be reached given that
# vaccination rate of 1% is administered.

Initial_values = c(S=989999, I=1, R=10000,V=10000,Total=1)
parameters = c(beta=0.6, sigma=0.25, v=10000)

# time
time2=seq(from=1, to=67, by=1)

# model equation-SIR vaccine model A ####

vm_model_A <-function(time2,state,parameters){
  with(as.list(c(state,parameters)),{
    N=S+I+R
    F_A = beta*(I/N)*(S/N)
    dS = -F_A*S - v
    dI = F_A*S - sigma*I
    dR = sigma*I + v
    dV = v
    dTotal = F_A*S 
    
    return(list(c(dS,dI,dR,dV,dTotal)))
  }
  )
}

#solving the equations
VM_A = as.data.frame(ode(y = Initial_values,
                         func = vm_model_A,
                         parms = parameters,
                         times = time2)) 
print(VM_A[67,])

Initial_values = c(S=328409.6, I=8.610055, R=671581.8, V=670000, Total=1590.414)
parameters = c(beta=0.6, sigma=0.25, v=0)

# time
time3=seq(from=67, to=500, by=1)

# model equation- SIR vaccine model A ####

vm2_model_A <-function(time3,state,parameters){
  with(as.list(c(state,parameters)),{
    N=S+I+R
    F_A = beta*(I/N)*(S/N)
    dS = -F_A*S - v
    dI = F_A*S - sigma*I
    dR = sigma*I + v
    dV = v
    dTotal = F_A*S 
    
    return(list(c(dS,dI,dR,dV,dTotal)))
  }
  )
}

#solving the equations
VM2_A = as.data.frame(ode(y = Initial_values,
                          func = vm2_model_A,
                          parms = parameters,
                          times = time3)) 

VMfinal_A = rbind(VM_A,VM2_A)
#View(VMfinal_A)
VMfinal_A = VMfinal_A[-68,]

# plot output SIR vaccine model A
N=1000000
VMfinal_A$Total = 100*VMfinal_A$Total/N
VMfinal_A$S = 100*VMfinal_A$S/N
VMfinal_A$I = 100*VMfinal_A$I/N
VMfinal_A$R = 100*VMfinal_A$R/N
VMfinal_A$V = 100*VMfinal_A$V/N

print(VMfinal_A[67,])
print(VMfinal_A[500,])

plot(VMfinal_A$time,
     VMfinal_A$Total,type = "l", col="blue",
     ylab = "Proportion, %", ylim = c(0,0.4),
     xlab = "Arbitrary time, t")
abline(v=67)
title("Vaccination Rate: 1% per Unit t", adj=0)
text(x=90, y=0.1, "HIT was achieved at t=67",cex=0.8, adj=0)

write.csv(VMfinal_A, file = "SIR_HIT_1%.csv")