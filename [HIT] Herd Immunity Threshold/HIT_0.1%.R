# This code aims to predict when HIT would be reached given that
# vaccination rate of 0.1% is administered.

Initial_values = c(S=998999, I=1, R=1000, V=1000, Total=1)
parameters = c(beta=0.6, sigma=0.=25, v=1000)

# time
time12=seq(from=1, to=134, by=1)

vm_model_A <-function(time12,state,parameters){
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
                         times = time12)) 
print(VM_A[134,])

Initial_values = c(S=417084.1 , I=0.3854066 , R=582915.5 ,V=1000 , Total=581915.9)
parameters = c(beta=0.6, sigma=0.25, v=0)

# time
time13=seq(from=134, to=200, by=1)

# model equation- SIR vaccine model A ####

vm2_model_A <-function(time13,state,parameters){
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
                          times = time13)) 

VMfinal_A = rbind(VM_A,VM2_A)
#View(VMfinal_A)
VMfinal_A = VMfinal_A[-135,]

# plot output SIR vaccine model A
N=1000000
VMfinal_A$Total = 100*VMfinal_A$Total/N
VMfinal_A$S = 100*VMfinal_A$S/N
VMfinal_A$I = 100*VMfinal_A$I/N
VMfinal_A$R = 100*VMfinal_A$R/N
VMfinal_A$V = 100*VMfinal_A$V/N

print(VMfinal_A[134,])
print(VMfinal_A[200,])

plot(VMfinal_A$time,
     VMfinal_A$Total,type = "l", col="blue",
     ylab = "Proportion, %", ylim = c(0,80),
     xlab = "Arbitrary time, t")
abline(v=134)
title("Vaccination Rate: 0.1% per Unit t", adj=0)
text(x=130, y=75, "HIT was achieved at t=134",cex=0.8, adj=1)


write.csv(VMfinal_A, file = "SIR_HIT_0.1%.csv")