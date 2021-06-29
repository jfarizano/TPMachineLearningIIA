library(caret)
library(rpart)
library(rpart.plot)
library(rattle)
setwd("~/TPMachineLearningIIA")

# Leemos los datos
yeast <- read.table("yeast.data", header = FALSE, col.names = c("sname", "mcg", "gvh", "alm", "mit", "erl", "pox", "vac", "nuc", "type"))
yeast <- yeast[, -1]

# Creamos los folds
indexData <- createFolds(t(yeast[, "type"]), k = 5)

# Creamos los conjuntos de test y entrenamiento
kElegido <- 3
yeastTest <- yeast[indexData[[kElegido]], ]
yeastTrain <- yeast[setdiff(seq(1:dim(yeast)[1]), indexData[[kElegido]]), ]

cpElegido <- 0.03

# -------------------------------------------------------------------------------------
# GANANCIA DE INFORMACIÓN
fitGain <- rpart(type~., data = yeastTrain,
              parms = list(split = 'information'))
# summary(fitGain)
fancyRpartPlot(fitGain, caption = NULL)
predictYeastGain <- predict(fitGain, yeastTest[,-9], type = 'class')
# confusionMatrix(predictYeast, yeastTest[,9])
matrixGain <- table(predictYeastGain, yeastTest[,9])
acuraccyGain <- sum(diag(matrixGain)) / sum(matrixGain)
precisionGain <- diag(matrixGain) / rowSums(matrixGain)
recallGain <- diag(matrixGain) / colSums(matrixGain)
measuresGain <- rbind(precisionGain, recallGain)

predictYeastGainTrain <- predict(fitGain, yeastTrain[,-9], type = 'class')
matrixGainTrain <- table(predictYeastGainTrain, yeastTrain[,9])
acuraccyGainTrain <- sum(diag(matrixGainTrain)) / sum(matrixGainTrain)

fitGainPruned <- prune(fitGain, cp=cpElegido)
# summary(fitGainPruned)
fancyRpartPlot(fitGainPruned, caption = NULL)
predictYeastGainPruned <- predict(fitGainPruned, yeastTest[,-9], type = 'class')
# confusionMatrix(predictYeastPruned, yeastTest[,9])
matrixGainPruned <- table(predictYeastGainPruned, yeastTest[,9])
acuraccyGainPruned <- sum(diag(matrixGainPruned)) / sum(matrixGainPruned)
precisionGainPruned <- diag(matrixGainPruned) / rowSums(matrixGainPruned)
recallGainPruned <- diag(matrixGainPruned) / colSums(matrixGainPruned)
measuresGainPruned <- rbind(precisionGainPruned, recallGainPruned)


# -------------------------------------------------------------------------------------
# GINI

fitGini <- rpart(type~., data = yeastTrain,
                 parms = list(split = 'gini'))
# summary(fitGini)
fancyRpartPlot(fitGini, caption = NULL)
predictYeastGini <- predict(fitGini, yeastTest[,-9], type = 'class')
# confusionMatrix(predictYeast, yeastTest[,9])
matrixGini <- table(predictYeastGini, yeastTest[,9])
acuraccyGini <- sum(diag(matrixGini)) / sum(matrixGini)
precisionGini <- diag(matrixGini) / rowSums(matrixGini)
recallGini <- diag(matrixGini) / colSums(matrixGini)
measuresGini <- rbind(precisionGini, recallGini)

predictYeastGiniTrain <- predict(fitGini, yeastTrain[,-9], type = 'class')
matrixGiniTrain <- table(predictYeastGiniTrain, yeastTrain[,9])
acuraccyGiniTrain <- sum(diag(matrixGiniTrain)) / sum(matrixGiniTrain)

fitGiniPruned <- prune(fitGini, cp=cpElegido)
# summary(fitGiniPruned)
fancyRpartPlot(fitGiniPruned, caption = NULL)
predictYeastGiniPruned <- predict(fitGiniPruned, yeastTest[,-9], type = 'class')
# confusionMatrix(predictYeastPruned, yeastTest[,9])
matrixGiniPruned <- table(predictYeastGiniPruned, yeastTest[,9])
acuraccyGiniPruned <- sum(diag(matrixGiniPruned)) / sum(matrixGiniPruned)
precisionGiniPruned <- diag(matrixGiniPruned) / rowSums(matrixGiniPruned)
recallGiniPruned <- diag(matrixGiniPruned) / colSums(matrixGiniPruned)
measuresGiniPruned <- rbind(precisionGiniPruned, recallGiniPruned)


# -------------------------------------------------------------------------------------
# Datos
acuraccy <- rbind(acuraccyGainTrain, acuraccyGain, acuraccyGainPruned, acuraccyGiniTrain, acuraccyGini, acuraccyGiniPruned)
precision <- rbind (precisionGain, precisionGainPruned, precisionGini, precisionGiniPruned)
recall <- rbind(recallGain, recallGainPruned, recallGini, recallGiniPruned)

print(acuraccy)
print(precision)
print(recall)