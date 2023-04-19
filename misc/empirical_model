##february 2023 modelling, with thickness (CRST paper)
library(bootstrap)
library(Metrics)
library(e1071)
library(car)
library(relaimpo)
library(caret)

wd <- "/data/"
setwd(wd)
ice <- read.csv("ice_dec17.csv")

#check normality of data:

#volume (square root) 1.084e-07
hist(ice$volume)
shapiro.test(ice$volume) #not normal
shapiro.test(sqrt(ice$volume)) #square root
shapiro.test(log(ice$volume)) #volume
shapiro.test((sign(ice$volume) * abs(ice$volume)^(1/3))) #cuberoot

#ram length - not a continuous distribution
hist(ice$ram_length)

#ram thickness - not a continuous distribution
hist(ice$ram_thickness)

#ram extent (log) 0.0021
hist(ice$ram_extent)
shapiro.test(ice$ram_extent) #not normal
shapiro.test(sqrt(ice$ram_extent)) #square root
shapiro.test(log(ice$ram_extent)) #log
shapiro.test((sign(ice$ram_extent) * abs(ice$ram_extent)^(1/3))) #cuberoot

#stress (cube root) 0.00007
hist(ice$sigma) 
shapiro.test(ice$sigma) #not normal
shapiro.test(sqrt(ice$sigma)) #square root
shapiro.test(log(ice$sigma)) #volume
shapiro.test((sign(ice$sigma) * abs(ice$sigma)^(1/3))) #cuberoot

#ram volume (log) >0.05
hist(ice$ram_vol) 
shapiro.test(ice$ram_vol) #not normal
shapiro.test(sqrt(ice$ram_vol)) #square root
shapiro.test(log(ice$ram_vol)) #log
shapiro.test((sign(ice$ram_vol) * abs(ice$ram_vol)^(1/3))) #cuberoot

#ram percent(log) > 0.05
hist(ice$ram_percent)
shapiro.test(ice$ram_percent) #not normal
shapiro.test(sqrt(ice$ram_percent)) #square root
shapiro.test(log(ice$ram_percent)) #log
shapiro.test((sign(ice$ram_percent) * abs(ice$ram_percent)^(1/3))) #cuberoot


##make a new df with the transformed variables
df <- matrix(ncol=7, nrow=nrow(ice))
df <- as.data.frame(df)
variables <- colnames(ice)[c(3,4,5,6,7,10,11)]
colnames(df) <- variables

df$volume <- sqrt(ice$volume) #square root
df$ram_length <- ice$ram_length
df$ram_thickness <- ice$ram_thickness
df$ram_extent <- log(ice$ram_extent)
df$sigma <- (sign(ice$sigma) * abs(ice$sigma)^(1/3)) #cuberoot
df$ram_vol <- log(ice$ram_vol)
df$ram_percent <- log(ice$ram_percent) #log

#get range of skewness in old versus new variables:
iterations <- length(variables)
vars = 3

output <- matrix(ncol=vars, nrow=iterations)

for (i in 1:length(variables)) {
  
  output[i,1] <- variables[i]
  output[i,2] <- round(skewness(ice[,variables[i]]),2)
  output[i,3] <- round(skewness(df[,variables[i]]),2)
  
}

colnames(output) <- c("variable", "before", "after")

#start the actual model

##only 2nd order interaction variables
full.model <- lm(sigma~(volume+ram_length+ram_thickness+ram_extent+
                           ram_vol+ram_percent)^2, data = df) 
summary(full.model)

#vif to test for multicollinearity
vif(lm(sigma ~ volume+ram_length+ram_thickness+ram_extent+ram_vol+
         ram_percent, data = df)) #remove ram_vol

vif(lm(sigma ~ volume+ram_length+ram_thickness+ram_extent+
         ram_percent, data = df)) #remove vol

vif(lm(sigma ~ ram_length+ram_thickness+ram_extent+
         ram_percent, data = df)) #all low VIF

##multicollinearity corrected model
full.model <- lm(sigma~(ram_length+ram_thickness+ram_extent+
                          ram_percent)^2, data = df) 
summary(full.model) #remove ram thickness :(

mod1 <- update(full.model, ~. - ram_thickness) 
summary(mod1)
AIC(full.model, mod1) 
anova(full.model, mod1) #remove ram_length:ram extent

mod2 <- update(mod1, ~. -ram_length:ram_extent) #back down
summary(mod2)
AIC(mod1, mod2)
anova(mod1, mod2) #remove ram_percent:ram_thickness

mod3 <- update(mod2, ~. - ram_percent:ram_thickness ) 
summary(mod3)
AIC(mod2, mod3)
anova(mod2, mod3) #remove ram_extent:ram_percent

mod4 <- update(mod3, ~. -ram_extent:ram_percent) 
summary(mod4)
AIC(mod3, mod4)   
anova(mod3, mod4) # remove ram_length:ram_thickness

mod5 <- update(mod4, ~. -ram_length:ram_thickness)
summary(mod5)
AIC(mod4, mod5)   
anova(mod4, mod5) #remove ram_percent

mod6 <- update(mod5, ~. -ram_percent) ##best model
summary(mod6)
AIC(mod5, mod6)   
anova(mod5, mod6)

#model ram thickness and length as categorical variables

##ANCOVA time

df$ram_length <- as.factor(df$ram_length)
df$ram_thickess <- as.factor(df$ram_thickness)

vif(lm(sigma ~ ram_length+ram_thickness+ram_extent+
         ram_percent, data = df))#all low VIF

##multicollinearity corrected model
full.model <- lm(sigma~(ram_length+ram_thickness+ram_extent+
                          ram_percent)^2, data = df) 
summary(full.model) #remove ram_length:ram_extent

mod1 <- update(full.model, ~. -ram_length:ram_extent) 
summary(mod1)
AIC(full.model, mod1) 
anova(full.model, mod1) #remove ram_thickness:ram_extent 

mod2 <- update(mod1, ~. -ram_thickness:ram_extent) #back down
summary(mod2)
AIC(mod1, mod2)
anova(mod1, mod2) # remove ram_thickness:ram_percent

mod3 <- update(mod2, ~. - ram_thickness:ram_percent) 
summary(mod3)
AIC(mod2, mod3)
anova(mod2, mod3) #remove ram_length:ram_percent

mod4 <- update(mod3, ~. -ram_length:ram_percent) 
summary(mod4)
AIC(mod3, mod4)   
anova(mod3, mod4) # remove ram_percent 

mod5 <- update(mod4, ~. -ram_percent)
summary(mod5)
AIC(mod4, mod5)   
anova(mod4, mod5) #remove ram_length:ram_thickness

mod6 <- update(mod5, ~. -ram_length:ram_thickness) ##best model
summary(mod6)
AIC(mod5, mod6)   
anova(mod5, mod6) #remove 


###no transforms whatsoever:
##only 2nd order interaction variables
full.model <- lm(sigma~(volume+ram_length+ram_thickness+ram_extent+
                          ram_vol+ram_percent)^2, data = df) 
summary(full.model)

#vif to test for multicollinearity
vif(lm(sigma ~ volume+ram_length+ram_thickness+ram_extent+ram_vol+
         ram_percent, data = df)) #remove ram_vol

vif(lm(sigma ~ volume+ram_length+ram_thickness+ram_extent+
         ram_percent, data = df)) #remove vol

vif(lm(sigma ~ ram_length+ram_thickness+ram_extent+
         ram_percent, data = df)) #all low VIF

##multicollinearity corrected model
full.model <- lm(sigma~(ram_length+ram_thickness+ram_extent+
                          ram_percent)^2, data = ice) 
summary(full.model) #remove ram extent:ram_percent :(

mod1 <- update(full.model, ~. -ram_extent:ram_percent) 
summary(mod1)
AIC(full.model, mod1) 
anova(full.model, mod1) #remove ram_thickness

mod2 <- update(mod1, ~. -ram_thickness) #back down
summary(mod2)
AIC(mod1, mod2)
anova(mod1, mod2) #remove ram_percent:ram_thickness

mod3 <- update(mod2, ~. - ram_percent:ram_thickness ) 
summary(mod3)
AIC(mod2, mod3)
anova(mod2, mod3) #remove ram_length:ram_extent

mod4 <- update(mod3, ~. -ram_length:ram_extent) 
summary(mod4)
AIC(mod3, mod4)   
anova(mod3, mod4) # remove ram_percent

mod5 <- update(mod4, ~. -ram_percent)
summary(mod5)
AIC(mod4, mod5)   
anova(mod4, mod5) #remove ram_extent:ram_thickness

mod6 <- update(mod5, ~. -ram_extent:ram_thickness) ##best model
summary(mod6)
AIC(mod5, mod6)   
anova(mod5, mod6)

#calculate relative importance of variables:
mod6Importance <- calc.relimp(mod6) 
mod6Importance@lmg #ram length imp, for example

###make a table for model statistics like this

model_stats <- as.data.frame(matrix(nrow = length(mod6$coefficients), ncol = 7))
model_stats[,1] <- c("intercept", attr(mod6$terms, "term.labels"))
model_stats[,2] <- round(mod6$coefficients,2)
model_stats[,3] <- round(mod6$coefficients,2)
model_stats[2:nrow(model_stats),3] <- mod6Importance@lmg 


round(model_stats$V3 * 100,2)


hist(ice$sigma)
abline(v = mean(ice$sigma), col = 'red')
abline(v= mean(ice$sigma) + 3 * sd(ice$sigma))










ci2d3 <- read.csv('/Users/jsmith/Desktop/Research/ci2d3_data.csv')
max(ci2d3$centroid_x)
min(ci2d3$centroid_x)
max(ci2d3$centroid_y)
min(ci2d3$centroid_y)

sqrt((max(ci2d3$centroid_y) - min(ci2d3$centroid_y))^2 +
(max(ci2d3$centroid_x) - min(ci2d3$centroid_x))^2 ) / 1000






mod1 <- lm(ice$sigma ~ ice$ram_length, data = ice)
summary(mod1)

xn <- (1:nrow(ice))


plot(ice$ram_length, ice$sigma)
abline(lm(ice$sigma ~ ice$ram_length), col = 'red')

plot(ice$ram_thickness, ice$sigma)
abline(lm(ice$sigma ~ ice$ram_thickness), col = 'red')


summary(lm(ice$sigma ~ ice$ram_thickness))








n <- length(df[, 1])
ind <- c(1:n)
RMSE <- NULL
MAE <- NULL
RSQ <- NULL 
COEF <- NULL

set.seed(123)
k <- 10
ksize <- n%/%k
data <- df
for (i in 1:(k - 1)) {
  ids <- c(1:length(data[, 1]))
  randids <- sample(ids, ksize)
  samp <- data[-randids, ]
  val <- data[randids, ]
  mod <- lm(sigma ~ ram_length+ram_extent+ ram_length:ram_percent+ram_extent:ram_thickness, data = samp)
  pred <- predict(mod, val)
  RMSE[i] <- rmse(pred, val$sigma)
  MAE[i] <- mae(pred, val$sigma)
  RSQ[i] <- summary(mod)$adj.r.squared
  #COEF[i] <- summary(mod)$coefficients
  data <- samp
}


hist(RMSE)
hist(MAE)
hist(COEF)
hist(RSQ)

mean(RMSE)
mean(MAE)
mean(COEF)
mean(RSQ)

train_control <- trainControl(method="cv", number=10) 
boot_model <- train(sigma ~ ram_length+ram_extent+ ram_length:ram_percent+ram_extent:ram_thickness, data = df,method = "lm", trControl=train_control)
print(boot_model)
