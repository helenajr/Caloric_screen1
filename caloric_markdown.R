## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- message= FALSE------------------------------------------------------------------------------------------------------------
library (tidyverse) #Load the tidyverse

data <- read_csv("caloric_tidy.csv")


## -------------------------------------------------------------------------------------------------------------------------------
data <- mutate(data, 
       CP = ((WR + CR)-(WL + CL))/(WR + WL + CR + CL)*100,
       DP = ((WR + CL) - (WL + CR))/(WR + WL + CR + CL)*100)


## -------------------------------------------------------------------------------------------------------------------------------
data <- mutate(data, 
                CP_call = if_else(((CP > -20) & (CP < 20)), "Neg", "Pos"),
                DP_call = if_else(((DP > -20) & (DP < 20)), "Neg", "Pos"),
                Hypo_call  = if_else(((WR < 8) & (WL < 8) & (CR < 8) & (CL <8)), "Hypo", "Normal"),
                gold_call = if_else(((CP_call == "Pos") | (DP_call == "Pos") | (Hypo_call == "Hypo")), 
                                    "Pos", "Neg"))



## -------------------------------------------------------------------------------------------------------------------------------
#Calculates monothermal screen score depending on order
data <- mutate(data,
               MS = if_else(Order == "W", (WR - WL)/(WR + WL)*100, (CR - CL)/(CR + CL)*100))

#Convert to non-negative number
data <- mutate(data,
               MS = if_else(MS <0, MS*-1, MS))

#Identifies patients where the monothermal screen has indicated bilateral hypofunction
data <- mutate(data,
                HS = if_else(((Order == "W") & (WR < 8) & (WL <8)), "Hypo",
                        if_else(((Order =="C") & (CR < 8) & (CL <8)), "Hypo", "Normal")))

#If HS == "Hypo" MS needs to be converted to NA (finding will be significant irrespective of score)
data <- mutate(data,
               MS = ifelse((HS == "Hypo"), NA, MS))



## -------------------------------------------------------------------------------------------------------------------------------
#Convert gold_call and Order to factors
data$gold_call <- as.factor(data$gold_call)
data$Order <- as.factor(data$Order)

#Remove rows with NAs from data - (34 rows)
data <- data[!(is.na(data$MS)),]


## -------------------------------------------------------------------------------------------------------------------------------
xtabs(~ gold_call + Order, data = data)


## -------------------------------------------------------------------------------------------------------------------------------
#Split MS into quartiles
q <- quantile(data$MS, probs=seq(0,1,0.25))

#Find numbers of pos and neg in each quartile, calculate proportion (probability) positive
table1 <- table(data$gold_call[data$MS< q[2]])
p1 <- table1[[2]]/(table1[[1]]+table1[[2]])

table2 <- table(data$gold_call[data$MS >=q[2] & data$MS < q[3]])
p2 <- table2[[2]]/(table2[[1]]+table2[[2]])

table3 <- table(data$gold_call[data$MS >= q[3] & data$MS < q[4]])
p3 <- table3[[2]]/(table3[[1]]+table3[[2]])

table4 <- table(data$gold_call[data$MS >= q[4] ])
p4 <- table4[[2]]/(table4[[1]]+table4[[2]])

#Store all of these in 1 object called probs
probs <- c(p1, p2, p3, p4)
#Turn the probabilities into log odds
logits <- log(probs/(1-probs))

#Caluclate median MS in each group

meds <- c( median(data$MS[ data$MS < q[2] ]),
           median(data$MS[ data$MS >=q[2] & data$MS < q[3] ]),
           median(data$MS[ data$MS >=q[3] & data$MS < q[4] ]),
           median(data$MS[ data$MS >= q[4] ]))

#Plot it
plot(meds, logits, xlab = "MS", ylab = "log-odds(pos|MS)", las = 1)


## -------------------------------------------------------------------------------------------------------------------------------

#Logistic regression using MS, Order and gold_call
logistic <- glm(gold_call ~ MS * Order, data = data, family = "binomial")
summary(logistic)


## -------------------------------------------------------------------------------------------------------------------------------
#Create a list of evenly spaced values between 0 and 170
x_MS <- seq(0, 170, 1)

#Save coefficients from the model
b0 <- logistic$coefficients[1]
MS <- logistic$coefficients[2]
OrderW <- logistic$coefficients[3]
Inter <- logistic$coefficients[4]

#Equations for group C
c_logits <- b0 + MS*x_MS
c_probs <- exp(c_logits)/(1 + exp(c_logits))

#Equations for group W
w_logits <- b0 + MS*x_MS + OrderW*1 + Inter*x_MS*1
w_probs <- exp(w_logits)/(1 + exp(w_logits))

#Plot it
plot(x_MS, c_probs, 
     ylim=c(0,1),
     type="l", 
     lwd=3, 
     lty=2, 
     col="turquoise2",
     xlab="MS", ylab="P(outcome)", main="Logistic functions for warm and cold groups")


# Add the line for people who are in the b group
lines(x_MS, w_probs, 
      type="l", 
      lwd=3, 
      lty=2, 
      col="orangered")



## -------------------------------------------------------------------------------------------------------------------------------
library(ROCit)

data$prediction <- logistic$fitted.values

logroc_results <- rocit(score = data$prediction, class = data$gold_call, negref = "Neg")

plot(logroc_results, values = F, YI = F, legend = F)

summary(logroc_results)
ciAUC(logroc_results)


## -------------------------------------------------------------------------------------------------------------------------------
gold_sum <- summary(as.factor(data$gold_call))

gold_sig <- gold_sum[names(gold_sum) == "Pos"]
gold_insig <- gold_sum[names(gold_sum) == "Neg"]
prevalence <- (gold_sig/(gold_sig + gold_insig))*100

cat("Prevalence of significant findings in study is ", prevalence, "%")



## -------------------------------------------------------------------------------------------------------------------------------

Positives <- prevalence
Negatives <- 100-prevalence

data2 <- tibble(Cutoff=logroc_results$Cutoff, 
                      Sensitivity=logroc_results$TPR, 
                      Specificity=1-logroc_results$FPR,
                FN=Positives - (Positives*Sensitivity),
                FP=Negatives - (Negatives*Specificity),
                TN=Negatives*Specificity,
                TP=Positives*Sensitivity,
                NPV=TN/(TN + FN),
                PPV=TP/(FP + TP),
                FT= TP + FP)

#Find rows with sensitivity closest to 0.9, 0.95, 0.98 and 0.99 and display in a table
list90 <- which(abs(data2$Sensitivity - 0.9)==min(abs(data2$Sensitivity - 0.9)))
row90 <- round(sum(list90)/length(list90))

list95 <- which(abs(data2$Sensitivity - 0.95)==min(abs(data2$Sensitivity - 0.95)))
row95 <- round(sum(list95)/length(list95))

list98 <- which(abs(data2$Sensitivity - 0.98)==min(abs(data2$Sensitivity - 0.98)))
row98 <- round(sum(list98)/length(list98))

list99 <- which(abs(data2$Sensitivity - 0.99)==min(abs(data2$Sensitivity - 0.99)))
row99 <- round(sum(list99)/length(list99))

data3 <- slice(data2, row90, row95, row98, row99)

library(knitr)

kable(data3, digits = 2)



## -------------------------------------------------------------------------------------------------------------------------------
#Find numbers of pos and neg in each bin - would be more efficient with a loop
table1 <- table(data$gold_call[data$prediction< 0.1])
p1 <- table1[[2]]/(table1[[1]]+table1[[2]])

table2 <- table(data$gold_call[data$prediction >= 0.1 & data$prediction < 0.2])
p2 <- table2[[2]]/(table2[[1]]+table2[[2]])

table3 <- table(data$gold_call[data$prediction >= 0.2 & data$prediction < 0.3])
p3 <- table3[[2]]/(table3[[1]]+table3[[2]])

table4 <- table(data$gold_call[data$prediction >= 0.3 & data$prediction < 0.4])
p4 <- table4[[2]]/(table4[[1]]+table4[[2]])

table5 <- table(data$gold_call[data$prediction >= 0.4 & data$prediction < 0.5])
p5 <- table5[[2]]/(table5[[1]]+table5[[2]])

table6 <- table(data$gold_call[data$prediction >= 0.5 & data$prediction < 0.6])
p6 <- table6[[2]]/(table6[[1]]+table6[[2]])

table7 <- table(data$gold_call[data$prediction >= 0.6 & data$prediction < 0.7])
p7 <- table7[[2]]/(table7[[1]]+table7[[2]])

table8 <- table(data$gold_call[data$prediction >= 0.7 & data$prediction < 0.8])
p8 <- table8[[2]]/(table8[[1]]+table8[[2]])

table9 <- table(data$gold_call[data$prediction >= 0.8 & data$prediction < 0.9])
p9 <- table9[[2]]/(table9[[1]]+table9[[2]])

table10 <- table(data$gold_call[data$prediction >= 0.9 ])
p10 <- table10[[2]]/(table10[[1]]+table10[[2]])

#Store all of these in 1 object called probs
cal_probs <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)

#Have the x values midway through each bin
cal_x <- seq(0.05, 0.95, 0.1)

#Plot it
plot(cal_x, cal_probs, xlab = "Model predictions", ylab = "Fraction of positives", las = 1)


lines(cal_x, cal_x, 
      type="l", 
      lwd=3, 
      lty=3, 
      col="grey")

lines(cal_x, cal_probs, 
      type="l", 
      lwd=2, 
      lty=1, 
      col="blue")


