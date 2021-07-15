library(tidyverse)
library(ROCit)
library(ggpubr)

data <- read_csv("caloric_tidy.csv") #Raw measurements provided by Vicki

#Calculate CP and DP
data <- mutate(data, 
               CP = ((WR + CR)-(WL + CL))/(WR + WL + CR + CL)*100,
               DP = ((WR + CL) - (WL + CR))/(WR + WL + CR + CL)*100)

#Calculate MS and HS, depending on which temperature was performed first
data <- mutate(data,
               MS = if_else(Order == "W", (WR - WL)/(WR + WL)*100, (CR - CL)/(CR + CL)*100),
               HS = if_else(((Order == "W") & (WR < 8) & (WL <8)), "Hypo",
                            if_else(((Order =="C") & (CR < 8) & (CL <8)), "Hypo", "Normal")))

#Convert to non-negative number - direction is not important here
data <- mutate(data,
               MS = if_else(MS <0, MS*-1, MS),
               CP = if_else(CP <0, CP*-1, CP),
               DP = if_else(DP <0, DP*-1, DP))

#Determine if there is a significant finding for the full test
data <- mutate(data, 
               CP_call = if_else((CP > 20), "Pos", "Neg"),
               DP_call = if_else((DP > 20), "Pos", "Neg"),
               Hypo_call  = if_else(((WR < 8) & (WL < 8) & (CR < 8) & (CL <8)), 
                                    "Hypo", "Normal"),
               gold_call = if_else((
                 (CP_call == "Pos") | (DP_call == "Pos") | (Hypo_call == "Hypo")), 
                 "Pos", "Neg"))

#If HS == "Hypo" MS is converted to Inf, this means it will count as being positive whatever MS cutoff is used
data <- mutate(data,
               MS = ifelse((HS == "Hypo"), Inf, MS))

#ROC analysis of both hot and cold together
roc_results <- rocit(score = data$MS, class = data$gold_call, negref = "Neg")

summary(roc_results)
ciAUC(roc_results)

plot(roc_results, values = F, YI = F)

data_tibble <- tibble(Cutoff=roc_results$Cutoff, 
                      TPR=roc_results$TPR, 
                      FPR=roc_results$FPR)


############ Cut off vs false positive rate

plot_data <- filter(data_tibble, Cutoff <=30)

plot1 <- ggplot(data = plot_data, aes(Cutoff, FPR)) + geom_point(size = 0.2)

plot_data <- mutate(plot_data,
                    FNR = 1-TPR)

plot2 <- ggplot(data = plot_data, aes(Cutoff, FNR)) + geom_point(size = 0.2)

figure0 <- ggarrange(plot1, plot2, labels = c("B1", "B2"))
figure0

###Sensitivity and specificity table

#Calculate prevalence of positives
gold_sum <- summary(as.factor(data$gold_call))

gold_sig <- gold_sum[names(gold_sum) == "Pos"]
gold_insig <- gold_sum[names(gold_sum) == "Neg"]
prevalence <- (gold_sig/(gold_sig + gold_insig))*100

cat("Prevalence of significant findings in study is ", prevalence, "%")

#Make a lovely table of the useful information
Positives <- prevalence
Negatives <- 100-prevalence

data2 <- tibble(Cutoff=roc_results$Cutoff, 
                Sensitivity=roc_results$TPR, 
                Specificity=1-roc_results$FPR,
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

data3 <- round(data3, 2)

write.table(data3, file = "summary1.txt", sep = ",", quote = FALSE, row.names = F)

##### COLD & WARM ####

c_data <- filter(data, Order == "C")
w_data <- filter(data, Order == "W")

croc_results <- rocit(score = c_data$MS, class = c_data$gold_call, negref = "Neg")
wroc_results <- rocit(score = w_data$MS, class = w_data$gold_call, negref = "Neg")

summary(croc_results)
ciAUC(croc_results)
summary(wroc_results)
ciAUC(wroc_results)

plot(croc_results, values = F, YIndex = F, col = c("blue", "grey"), legend = F)
lines(wroc_results$TPR~wroc_results$FPR, col = "red", lwd = 1.75)
legend("bottomright", c("Warm test","Cool test", "Chance line"),
                                lty = c(1,1,2), col = c("red", "blue", "grey"))

##Making a plot with confidence intervals - doesn't work, I think because of the infinity values
# ci_warm <- ciROC(wroc_results, level = 0.95)
# ci_cold <- ciROC(croc_results, level = 0.95)
# 
# plot(ci_warm, legend = F, col = 2)
# lines(ci_cold$TPR~ci_cold$FPR, col = "blue", lwd = 1.5)
# lines(ci_cold$LowerTPR~ci_cold$FPR, col = "blue", lty = 2)
# lines(ci_cold$UpperTPR~ci_cold$FPR, col = "blue", lty = 2)
# legend("bottomright", c("Warm test",
#                         "Warm test 95% CI",
#                         "Cool test",
#                         "Cool test 95% CI"),
#        lty = c(1,2,1,2), col = c(2, 2, "blue", "blue"))


c_data_tibble <- tibble(Cutoff=croc_results$Cutoff, 
                      TPR=croc_results$TPR, 
                      FPR=croc_results$FPR,
                      Temp = "Cold")

w_data_tibble <- tibble(Cutoff=wroc_results$Cutoff, 
                        TPR=wroc_results$TPR, 
                        FPR=wroc_results$FPR,
                        Temp = "Warm")

## Making the summary tables
data2 <- tibble(Cutoff=wroc_results$Cutoff, 
                Sensitivity=wroc_results$TPR, 
                Specificity=1-wroc_results$FPR,
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

data3 <- round(data3, 2)

write.table(data3, file = "summary_warm.txt", sep = ",", quote = FALSE, row.names = F)

data2 <- tibble(Cutoff=croc_results$Cutoff, 
                Sensitivity=croc_results$TPR, 
                Specificity=1-croc_results$FPR,
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

data3 <- round(data3, 2)

write.table(data3, file = "summary_cold.txt", sep = ",", quote = FALSE, row.names = F)



####### Putting both cold and warm on the same plot #####

combo <- bind_rows(c_data_tibble, w_data_tibble)


plot_combo <- filter(combo, Cutoff <=30)

p1 <- ggplot(data = plot_combo, aes(Cutoff, FPR)) + 
  geom_point(size = 0.2, aes(color = Temp))

p1 <- p1 + scale_color_manual(values = c("#00BFC4", "#F8766D"))

plot_combo <- mutate(plot_combo,
                    FNR = 1-TPR)

p2 <- ggplot(data = plot_combo, aes(Cutoff, FNR)) + geom_point(size = 0.2, aes(color = Temp))

p2 <- p2 + scale_color_manual(values = c("#00BFC4", "#F8766D"))

figure <- ggarrange(p1, p2, labels = c("B1", "B2"), common.legend = TRUE)
figure

