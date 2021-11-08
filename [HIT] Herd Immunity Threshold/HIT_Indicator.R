# This code aims to show how HIT predicts the total infections for the entire 
# duration of epidemic given that no vaccination is administered.

library(deSolve)
dev.off()

Initial_values = c(S=999999, I=1, R=0, Total=1)
parameters = c(beta=0.6, sigma=0.25)
  
# time
time=seq(from=1, to=150, by=1)

sir_model_A <-function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    N=S+I+R
    F_A = beta*(I/N)*(S/N)
    dS = -F_A*S
    dI = F_A*S - sigma*I
    dR = sigma*I
    dTotal = F_A*S 
    
    return(list(c(dS,dI,dR,dTotal)))
  }
  )
}

#solving the equations
SIR_model_A = as.data.frame(ode(y =Initial_values,
                                func = sir_model_A,
                                parms = parameters,
                                times = time1)) 

# plot output SIR model 
N=1000000
SIR_model_A$S = 100*SIR_model_A$S/N
SIR_model_A$I = 100*SIR_model_A$I/N
SIR_model_A$R = 100*SIR_model_A$R/N
SIR_model_A$Total = 100*SIR_model_A$Total/N

print(SIR_model_A[195,])
print(SIR_model_A[200,])

plot(SIR_model_A$time,
     SIR_model_A$S,type = "l", col="blue",
     ylab = "Proportion, %", ylim = c(0,100),
     xlab = "Arbitrary time, t")
lines(SIR_model_A$I, col="red")
lines(SIR_model_A$R, col="green")
lines(SIR_model_A$Total, col="Black")
abline(h=58.33,lty = 2)
text(x=25, y=90, "S(t)",cex=0.8)
text(x=50, y=15, "I(t)",cex=0.8)
text(x=50, y=28, "R(t)",cex=0.8)
text(x=120, y=63, "Total Infections",cex=0.8)
text(x=5, y=55, "HIT",cex=0.8)

write.csv(SIR_model_A, file = "SIR_HIT_Indicator.csv")
#*****************************************************