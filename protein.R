# protein.R
# classification of a protein data set
# D.McNeall

# ---------------------------------------------------------
# 0. Packages, palettes and setup
# ---------------------------------------------------------
set.seed(42)

library(RColorBrewer)
library(e1071)
library(MASS)

paired <- brewer.pal(6,'Paired') 
blues <- brewer.pal(9,'Blues')

mitoblue <- '#A2D9F7'
pogreen  <- '#A8D4AF'
ersalmon <- '#F6B096'

catpat <- paired
catpat[1] <- ersalmon
catpat[2] <- mitoblue
catpat[6] <- pogreen
catpat[3] <- paired[6]
catpat[4] <- paired[2]
catpat[5] <- paired[4]

# ---------------------------------------------------------
# 1. data
# ---------------------------------------------------------
protein <- read.csv('protein.csv')
# Or, we could just load up the truncated file (added 10th December 2015)
# This breaks the colour palettes at the moment.
protein_trunc <- read.csv('protein_trunc.csv')

# factors are alphabetical, so as.numeric gives
# ER          = 1
# MIT0        = 2
# MITO/ER     = 3
# MITO/PO     = 4
# MITO/PO/ER  = 5
# PO          = 6

# Training data 
train <- subset(protein, select = c('tmd_gravy', 'tail_charge'))
cl <- protein$location

# Data truncated to just the main 3 classes


#protein_trunc <- subset(protein, location=='ER' | location=='MITO' | location=='PO')

train_trunc <- subset(protein_trunc, select = c('tmd_gravy', 'tail_charge'))
cl_trunc <- protein_trunc$location

upper <- apply(train_trunc, 2, max)
lower <- apply(train_trunc, 2, min)

tmd_gravy_seq <- seq(from=lower[1], to=upper[1], length.out = 50)
tail_charge_seq <- seq(from=lower[2], to=upper[2], length.out = 50)

# A grid of test points that covers the entire data space
test <- as.matrix(expand.grid(tmd_gravy_seq, tail_charge_seq))
colnames(test) <- colnames(train_trunc)

# skip the last row, as there is no tail charge 
target<- read.table('protein_target.txt', sep = '\t', nrows = 244, header = TRUE)
colnames(target) <- colnames(train_trunc)

mutants <- read.csv('mutants.csv')

# A little aside to get the factors in the data sets to be the same. rbind()
# changes the as.numeric factors to match. 
combined <- rbind(protein, mutants)

mutants.ix <- seq(from = (nrow(protein)+1), to = nrow(protein)+nrow(mutants), by = 1)
mutants_refactor <- combined[mutants.ix, ]

# "Shared location" data
protein_shared <- subset(protein, location=='MITO/ER' | location=='MITO/PO' | location=='MITO/PO/ER')


# some new test sites

test.sites <- data.frame(cbind(c(1.54, 1.43, 2.34), c(3.9, 4.6, -0.1)))

colnames(test.sites) <- c('tmd_gravy', 'tail_charge')

#1) TMD GRAVY = 1.54 Charge =  3.9
#2) TMD GRAVY = 1.43 Charge =  4.6
#3) TMD GRAVY = 2.34 Charge = -0.1


# ---------------------------------------------------------
# Initial plot of the data
# ---------------------------------------------------------

dev.new()
par(las = 1)
plot(protein$tmd_gravy, protein$tail_charge, col = catpat[protein$location], pch = 19)
points(mutants_refactor$tmd_gravy, mutants_refactor$tail_charge, col = catpat[mutants_refactor$location], pch = 3)

legend('topleft', legend = unique(protein$location), col = catpat[unique(protein$location)], pch = 20, text.col = catpat[unique(protein$location)] )
dev.print(pdf, file = 'protein_all.pdf', width = 6, height = 6)

# ---------------------------------------------------------
# All the data and the mutant data together
# ---------------------------------------------------------
dev.new()
par(las = 1)
plot(protein$tmd_gravy, protein$tail_charge, bg = catpat[protein$location],
 pch = 21, col = 'black',
 xlab = 'tmd_gravy', ylab = 'tail_charge')

points(mutants_refactor$tmd_gravy, mutants_refactor$tail_charge, col = 'black', bg = catpat[mutants_refactor$location], pch = 24)

legend('topleft', legend = unique(protein$location), pt.bg = catpat[unique(protein$location)], pch = 21, text.col = catpat[unique(protein$location)], cex = 0.8, bty = 'n' )
legend('top', legend = c( 'triangles indicate mutants'), pch = c(24), cex = 0.8, bty = 'n')
dev.print(pdf, file='protein_and_mutants.pdf', width = 6, height = 6)




# Data truncated to just the 3 main classes
dev.new()
par(las = 1, fg = 'grey')
plot(protein_trunc$tmd_gravy, protein_trunc$tail_charge, col = catpat[protein_trunc$location], pch = 19)

legend('topleft', legend = unique(protein_trunc$location), col = catpat[unique(protein_trunc$location)], pch = 20, text.col = catpat[unique(protein_trunc$location)] )
dev.print(pdf, file='protein_trunc.pdf', width = 6, height = 6)

# ---------------------------------------------------------
# Support vector machine classification
# ---------------------------------------------------------

# Non-probabilistic SVM fit
svm.fit <- svm(location~., data = protein_trunc, probability = FALSE)
svm.pred <- predict(svm.fit, newdata=test, probability = FALSE)

# misclassification rate for the svm is ~18%
svm.mcr = 1 - (sum(predict(svm.fit)==protein_trunc$location)/ length(protein_trunc$location))

dev.new()
par(las = 1)
plot(test, col = catpat[svm.pred], pch = 20, cex = 0.6, main = 'SVM classifier prediction')
points(protein_trunc$tmd_gravy, protein_trunc$tail_charge, col = 'black',bg = catpat[protein_trunc$location], pch = 21)

legend('topleft', legend = unique(protein_trunc$location),
col = catpat[unique(protein_trunc$location)],
pch = 20, 
text.col=catpat[unique(protein_trunc$location)],
 bg = 'white')
dev.print(pdf, file = 'svm_prediction.pdf', width = 6, height = 6)
 
 
# --------------------------------------------------------- 
# Probabilistic SVM fit 
# ---------------------------------------------------------
# Why does this predict some non-seen classes?
svm.fit.prob <- svm(location~., data = protein_trunc, probability = TRUE)
svm.pred.prob <- predict(svm.fit.prob, newdata=test, probability = TRUE)

# How many does the svm fit correctly?

fit.hit <- sum(svm.fit.prob$fitted == protein_trunc$location)
fit.miss <- length(protein_trunc$location) - fit.hit
prob.mcr <- fit.miss / length(protein_trunc$location) 

svm.prob.df <- attr(svm.pred.prob,'prob')
mitomat <- matrix(svm.prob.df[,'MITO'], nrow = length(tmd_gravy_seq))
pomat <- matrix(svm.prob.df[,'PO'], nrow = length(tmd_gravy_seq))
ermat <- matrix(svm.prob.df[,'ER'], nrow = length(tmd_gravy_seq))


protein_trunc.fit <- cbind(protein_trunc,as.character(svm.fit.prob$fitted))
colnames(protein_trunc.fit)[4] <- c('predicted')
write.matrix(protein_trunc.fit, file = 'protein_trunc_fit.txt')



n <- length(protein_trunc$location)
cv.out <- rep(NA, n)
# A simple looping leave-one-out cross validation
for(i in 1:n){
	
	cv.fit <-  svm(location~., data = protein_trunc[-i, ], probability = TRUE)
	cv.pred <- predict(cv.fit, newdata=protein_trunc[i, ], probability = TRUE)
	cv.out[i] <- as.character(cv.pred) 

}

loo.hitrate<- sum(cv.out == as.character(protein_trunc$location)) / length(protein_trunc$location)

labcex = 0.6
pt.cex = 1.2
dev.new()
par(las=1, mar = c(4,4,1,1))
plot(test, col = catpat[svm.pred], type = 'n', main = '',
xlab = 'TMD GRAVY', ylab = 'Tail charge')#, main = 'SVM classifier class probability')

contour(tmd_gravy_seq, tail_charge_seq, mitomat,
	add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[2],labcex = labcex)
contour(tmd_gravy_seq, tail_charge_seq, pomat,
	add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[6], labcex = labcex)
contour(tmd_gravy_seq, tail_charge_seq, ermat,
	add = TRUE,levels = c(0.5,0.6,0.7,0.8), col = catpat[1], labcex = labcex)
    
points(protein_trunc$tmd_gravy, protein_trunc$tail_charge,
	col = 'black',bg = catpat[protein_trunc$location], pch = 21, cex = pt.cex)
points(test.sites$tmd_gravy, test.sites$tail_charge,
	col = 'black', pch = 21, lwd = 1.5, cex = pt.cex)

legend('topleft', legend = c(as.character(unique(protein_trunc$location)), 'test'),
	col = 'black', pt.bg = c(catpat[unique(protein_trunc$location)],'white'), pch = 21, text.col = c(catpat[unique(protein_trunc$location)], 'black'), bg = 'white', pt.lwd = c(rep(1,3),1.5), pt.cex = pt.cex)
dev.print(pdf, file='svm_probability_contour.pdf', width = 5, height = 5)

dev.print(tiff, file='svm_probability_contour.tif', width = 480, height = 480)



svm.protein <- svm(location~., data = protein, probability = TRUE)#!! 
protein.hit <- sum(svm.protein$fitted == protein$location)

protein.fit <- cbind(protein,as.character(svm.protein$fitted))
colnames(protein.fit)[4] <- c('predicted')

write.matrix(protein.fit, file = 'protein_fit.txt')


# ---------------------------------------------------------
# Predict the location in the target data
# ---------------------------------------------------------

target.prob <- predict(svm.fit.prob, newdata = target, probability = TRUE )

dev.new()
par(las=1)
plot(target, type = 'n',  main = 'SVM classifier target predictions')
points(target, col = catpat[target.prob], cex = 0.8, pch = 21)
contour(tmd_gravy_seq, tail_charge_seq, mitomat, add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[2])
contour(tmd_gravy_seq, tail_charge_seq, pomat,add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[6])
contour(tmd_gravy_seq, tail_charge_seq, ermat,add = TRUE, levels = c(0.5,0.6,0.7,0.8), col = catpat[1])
points(protein_trunc$tmd_gravy, protein_trunc$tail_charge, col = 'black',bg = catpat[protein_trunc$location], pch = 21)
legend('topright', legend = unique(protein_trunc$location), col = catpat[unique(protein_trunc$location)], pch = 20, text.col = catpat[unique(protein_trunc$location)], bg = 'white')
dev.print(pdf, file = 'protein_target.pdf', width = 6, height = 6)


write.matrix(cbind(target,attr(target.prob, 'probabilities'),target.prob), file = 'target_predictions.txt')

# A diagram which has the "mutants" data on top of the 
# probability contours for the three unique classes.
dev.new(width = 6, height = 6)
par(las=1)
plot(target, type = 'n',  main = '', xlim = c(0.95,2.8), ylim = c(-3.2,5.6))

contour(tmd_gravy_seq, tail_charge_seq, mitomat, add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[2])
contour(tmd_gravy_seq, tail_charge_seq, pomat,add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[6])
contour(tmd_gravy_seq, tail_charge_seq, ermat,add = TRUE, levels = c(0.5,0.6,0.7,0.8), col = catpat[1])


points(protein_trunc$tmd_gravy, protein_trunc$tail_charge, col =  catpat[protein_trunc$location], pch = 21)
points(mutants_refactor$tmd_gravy, mutants_refactor$tail_charge, col = 'black', bg = catpat[mutants_refactor$location], pch = 24)

text(mutants_refactor$tmd_gravy, mutants_refactor$tail_charge, col =  catpat[mutants_refactor$location], labels = letters[1:length(mutants_refactor$location)], pos = 2)

legend('bottomright', legend = unique(protein$location), text.col = catpat[unique(protein$location)],cex = 0.8, text.font = 2, bty = 'n')
par(xpd = TRUE)
#legend('topright', legend = unique(protein$location), text.col = catpat[unique(protein$location)])
legend('top', legend = 'denotes mutant', text.col = 'black', pch = 24, pt.bg = c( 'grey'), inset = c(0, -0.1), bty = 'n')

dev.print(pdf, file='mutants_foreground.pdf', width = 6, height = 6)

# A diagram which has the "shared location" data on top of the 
# probability contours for the three unique classes.
dev.new(width = 6, height = 6)
par(las=1)
plot(target, type = 'n',  main = '', xlim = c(0.95,2.8), ylim = c(-3.2,5.6))

contour(tmd_gravy_seq, tail_charge_seq, mitomat, add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[2])
contour(tmd_gravy_seq, tail_charge_seq, pomat,add = TRUE, levels = c(0.5,0.6,0.7,0.8,1), col = catpat[6])
contour(tmd_gravy_seq, tail_charge_seq, ermat,add = TRUE, levels = c(0.5,0.6,0.7,0.8), col = catpat[1])


points(protein_trunc$tmd_gravy, protein_trunc$tail_charge, col =  catpat[protein_trunc$location], pch = 21)
points(protein_shared$tmd_gravy, protein_shared$tail_charge, col = 'black', bg = catpat[protein_shared$location], pch = 22)

text(protein_shared$tmd_gravy, protein_shared$tail_charge, col =  catpat[protein_shared$location], labels = letters[1:length(protein_shared$location)], pos = 2)

legend('bottomright', legend = unique(protein$location), text.col = catpat[unique(protein$location)],cex = 0.8, text.font = 2, bty = 'n')
par(xpd = TRUE)
#legend('topright', legend = unique(protein$location), text.col = catpat[unique(protein$location)])
legend('top', legend = 'denotes shared location', text.col = 'black', pch = 22, pt.bg = c( 'grey'), inset = c(0, -0.1), bty = 'n')

dev.print(pdf, file='shared_foreground.pdf', width = 6, height = 6)

stop()

dev.new()
filled.contour(tmd_gravy_seq, tail_charge_seq, mitomat,
 levels=seq(from=0, to=1, by=0.2) ,
 xlab='tmd_gravy', ylab='tail_charge', main='SVM MITO probability', col=blues)
 
dev.new()
filled.contour(tmd_gravy_seq, tail_charge_seq, pomat,
 levels = seq(from = 0, to = 1, by = 0.2) ,
 xlab='tmd_gravy', ylab='tail_charge', main='SVM PO probability', col=blues)

dev.new()
filled.contour(tmd_gravy_seq, tail_charge_seq, ermat,
 levels = seq(from = 0, to = 1, by = 0.2) ,
 xlab='tmd_gravy', ylab='tail_charge', main='SVM ER probability', col=blues)


# K-nearest-neighbours test classification.


knn.cv.summ <- function(train, cl, num.k){
	# function loop through k to get a quick summary of knn
	# misclassification rate
	misclass.no.vec <- rep(NA, num.k)
	
	for(i in 1:num.k){	
		cv <- knn.cv(train = train,cl = cl, k = i, prob = FALSE)
		result <- cv == cl
		misclass.no.vec[i] <- length(cl) - sum(result)
		misclass.rate.vec <- misclass.no.vec / length(cl)			
	}
misclass.rate.vec	
}

knn.summary <- knn.cv.summ(train = train, cl = cl, num.k = 15) 
knn.summary.trunc <- knn.cv.summ(train = train_trunc, cl = cl_trunc, num.k = 15) 




pdf(file = 'knn_misclass.pdf', width = 6, height = 6)
par(las = 1)
plot(knn.summary, type = 'o', main = 'knn misclassification rate',
ylim = c(0,1), xlab = 'nearest neighbours', ylab = 'misclass rate')
lines(knn.summary.trunc, type = 'o', col = 'red')
dev.off()

# So, for the truncated data set, k = 3 looks like a good bet.



# Non-probabilistic class prediction to start with
pred <- knn(train=train_trunc, cl=cl_trunc, test=test, k=3)

pdf(file = 'protein_trunc_knn.pdf', width = 6, height = 6)
par(las = 1)
plot(test, col = catpat[pred], pch = 20, cex = 0.6, main = 'K-nearest-neighbour classifier (k = 3)')
points(protein_trunc$tmd_gravy, protein_trunc$tail_charge, col = 'black',bg = catpat[protein_trunc$location], pch = 21)

legend('topleft', legend = unique(protein_trunc$location), col = catpat[unique(protein_trunc$location)], pch = 20, text.col = catpat[unique(protein_trunc$location)], bg = 'white')
dev.off()

# probabilistic knn
# WARNING, THIS SECTION IS INCORRECT AS CODED (but might be useful)
pred.prob <- knn(train=train_trunc, cl=cl_trunc, test=test, k=3, prob = TRUE)

MITO.ix <- pred.prob=='MITO'
PO.ix   <- pred.prob=='PO'
ER.ix   <- pred.prob=='ER'

anti <- 1-attr(pred.prob, 'prob')

mitoprob <- anti
poprob   <- anti
erprob   <- anti

mitoprob[MITO.ix] <- attr(pred.prob, 'prob')[MITO.ix]
poprob[PO.ix] <- attr(pred.prob, 'prob')[PO.ix]
erprob[ER.ix] <- attr(pred.prob, 'prob')[ER.ix]

mitomat <- matrix(mitoprob, nrow = length(tmd_gravy_seq))
pomat   <- matrix(poprob, nrow = length(tmd_gravy_seq))
ermat   <- matrix(erprob, nrow = length(tmd_gravy_seq))


par(las = 1)
plot(test, col = catpat[pred], pch = 20, cex = 0.6, main = 'K-nearest-neighbour classifier (k = 3)')
points(protein_trunc$tmd_gravy, protein_trunc$tail_charge, col = 'black',bg = catpat[protein_trunc$location], pch = 21)

contour(tmd_gravy_seq, tail_charge_seq, mitomat, add = TRUE, levels = c(0, 0.5,1))

legend('topleft', legend = unique(protein_trunc$location), col = catpat[unique(protein_trunc$location)], pch = 20, text.col = catpat[unique(protein_trunc$location)], bg = 'white')

filled.contour(tmd_gravy_seq, tail_charge_seq, mitomat, levels = seq(from = 0, to = 1, by = 0.2) , col = blues)

levels = c(0,0.334, 0.667, 1)





