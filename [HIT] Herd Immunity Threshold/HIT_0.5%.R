# This code aims to predict when HIT would be reached given that
# vaccination rate of 0.5% is administered.

Initial_values = c(S=994999, I=1, R=5000, V=5000, Total=1)
parameters = c(beta=0.6, sigma=0.25, v=5000)

# time
time6=seq(from=1, to=131, by=1)

vm_model_A <-function(time6,state,parameters){
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
                         times = time6)) 
print(VM_A[131,])

Initial_values = c(S=202869.8, I=0.2832988, R=797130, V=655000, Total=142130.2)
parameters = c(beta=0.6, sigma=0.25, v=0)

# time
time7=seq(from=131, to=300, by=1)

# model equation- SIR model A ####

vm2_model_A <-function(time7,state,parameters){
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
                          times = time7)) 

VMfinal_A = rbind(VM_A,VM2_A)
#View(VMfinal_A)
VMfinal_A = VMfinal_A[-132,]

# plot output SIR vaccine model A
N=1000000
VMfinal_A$Total = 100*VMfinal_A$Total/N
VMfinal_A$S = 100*VMfinal_A$S/N
VMfinal_A$I = 100*VMfinal_A$I/N
VMfinal_A$R = 100*VMfinal_A$R/N
VMfinal_A$V = 100*VMfinal_A$V/N

print(VMfinal_A[131,])
print(VMfinal_A[300,])

plot(VMfinal_A$time,
     VMfinal_A$Total,type = "l", col="blue",
     ylab = "Proportion, %", ylim = c(0,20),
     xlab = "Arbitrary time, t")
abline(v=131)
title("Vaccination Rate: 0.5% per Unit t", adj=0)
text(x=135, y=5.0, "HIT was achieved at t=131",cex=0.8, adj=0)


write.csv(VMfinal_A, file = "SIR_HIT_0.5%.csv")