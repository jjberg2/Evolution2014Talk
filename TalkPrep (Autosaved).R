setwd ( "/Users/JeremyBerg/Documents/Academics/MultiTraitPolygenicAdaptation/FinnishSoftware/driftsel/R")
sapply(list.files(pattern="*.R"), source, .GlobalEnv)
setwd("../data")
sapply(list.files(pattern="*.rda"), load, .GlobalEnv)



# genotype data
setwd("~/Documents/Academics/MultiTraitPolygenicAdaptation/PopData/Height")
temp <- read.delim ( "height.snp.genos",h=F,sep = " " , stringsAsFactors=F)
haploid.genos <- temp[,6:ncol(temp)]
diploid.genos <- matrix ( NA , nrow = 157 , ncol = 2838 )
#diploid.genos [ , 1:5 ] <- as.matrix(temp [ , 1:5])
for ( i in 1:2838 ) {
	diploid.genos [ , i ] <- haploid.genos [ , (i-1)*2+1 ] + haploid.genos [ , (i-1)*2+2 ]
	
	
}
my.data <- cbind ( temp[,1:5],diploid.genos)
diploid.genos.height <- my.data
new.diploid.genos.height<-diploid.genos.height[!duplicated ( diploid.genos.height$V2),]
diploid.genos.height<-new.diploid.genos.height
save ( diploid.genos.height , file = "height.diploid.genos")
diploid.genos.height <- diploid.genos.height [ order ( diploid.genos.height$V2 ) ,  ]

# load genotype data
load ("height.diploid.genos")



# pca data
.pcadir <- "/Users/JeremyBerg/Documents/Academics/MultiTraitPolygenicAdaptation/PopData/"
pcas<-read.table(paste(.pcadir,"euro_nooutlier.pcavec",sep=""),skip=1,as.is=TRUE)
eig.vals<-read.table(paste(.pcadir,"euro_nooutlier.pcaval",sep=""),skip=1,as.is=TRUE)



# gwas.data
height.gwas <- read.delim ( "europe.height.txt" )
height.gwas <- height.gwas[height.gwas$SNP %in% diploid.genos.height$V2,]
height.gwas <- height.gwas [ order ( height.gwas$SNP) , ]


# figure out what's going on with alleles UGH!
check.vec <- character ( nrow ( height.gwas ) )
for ( i in 1 : nrow ( height.gwas ) ) {
	
	if ( (height.gwas$A1 [ i ] == diploid.genos.height$V4 [ i ]) & (height.gwas$A2 [ i ] == diploid.genos.height$V5 [ i ]) ){
		check.vec [ i ] <- T
	} else if ( (height.gwas$A1 [ i ] == diploid.genos.height$V5 [ i ]) & (height.gwas$A2 [ i ] == diploid.genos.height$V4 [ i ]) ) {
		check.vec [ i ] <- F
	} else {
		check.vec [ i ] <- "STRAND!!!!!"
	}
}

#deal with a strand problem for one SNP
diploid.genos.height [ 109 , 4 ] <- "T"
diploid.genos.height [ 109 , 5 ] <- "C"






## flip alleles
for ( i in 1 : nrow ( height.gwas ) ) {
	
	if ( (height.gwas$A1 [ i ] == diploid.genos.height$V4 [ i ]) & (height.gwas$A2 [ i ] == diploid.genos.height$V5 [ i ]) ){
		next
	} else if ( (height.gwas$A1 [ i ] == diploid.genos.height$V5 [ i ]) & (height.gwas$A2 [ i ] == diploid.genos.height$V4 [ i ]) ) {
		temp1 <- diploid.genos.height$V4 [ i ]
		temp2 <- diploid.genos.height$V5 [ i ]
		diploid.genos.height$V4 [ i ] <- temp2
		diploid.genos.height$V5 [ i ] <- temp1
		diploid.genos.height [ i , 6 : ncol ( diploid.genos.height ) ] <- 2 - diploid.genos.height [ i , 6 : ncol ( diploid.genos.height ) ]
	}
}
save ( diploid.genos.height , file = "height.diploid.genos.als.flipped")
load ( "height.diploid.genos.als.flipped")



# .fam file
fam.file <- read.delim ( "../POPRES.inds" , h=F , sep = " ")
inds.to.keep<-as.numeric(unlist ( lapply ( strsplit ( pcas[,1] , ":") , function ( x ) x [1 ])))
indexes.to.keep <- 5+which ( fam.file [ , 1 ] %in% inds.to.keep )
reduced.dip.gens <- diploid.genos.height [ , c ( 1:5 , indexes.to.keep ) ]
mean.freqs <- rowMeans ( reduced.dip.gens [ , 6:ncol ( reduced.dip.gens ) ] ) / 2
var <- sum(mean.freqs*(1-mean.freqs)*height.gwas$EFF^2)

POPRES.gen.heights <- colSums(reduced.dip.gens[ , 6:ncol ( reduced.dip.gens ) ] * height.gwas$EFF)

my.vars <- numeric ( 10)
for ( i in 1:10 ) {
my.vars [ i ] <- (sum(POPRES.gen.heights*pcas[,i+1]))^2/(var*eig.vals[i,1])
}

my.heights <- POPRES.gen.heights - mean ( POPRES.gen.heights )

my.vars <- numeric ( 10)
for ( i in 1:10 ) {
my.vars [ i ] <- (sum(my.heights*pcas[,i+1]))^2/(var*eig.vals[i,1])
}

ind.vars <- ( t ( pcas [ , 2:11 ] ) %*% my.heights )^2
my.vars <- numeric ( 10 )
for ( i in 1:10 ) {
	my.vars [ i ] <- ind.vars [ i , 1 ] / sum ( ind.vars [ , 1 ] ) / ( eig.vals [ i , 1 ] / sum ( eig.vals [ 1:10 , 1 ] ) )
}



ind.vars <- ( t ( pcas [ , 2:11 ] ) %*% my.heights )^2
my.vars <- numeric ( 10 )
for ( i in 1:10 ) {
	my.vars [ i ] <- ind.vars [ i , 1 ] / ( sum ( ind.vars [ -i , 1 ] ) * sum(eig.vals[1:10,1])/sum(eig.vals[-i,1]) )
}


######
plot ( NA , xlim = c ( 0 , 60 ) , ylim = c ( 0 , 0.5 ) , bty = "n" , xlab = expression ( Q[PC1:2]/F[PC1:2] ) , ylab = "Density" )
curve ( dchisq ( x , 2 ) , col = "black" , add = T , lwd = 2 )
arrows ( sum ( my.vars [1:2]) , 0.5 , sum ( my.vars [1:2]) , 0 , col = 'red' , lwd = 2 )


pdf ( "../../../Evolution2014/Figs/HeightPCs12.pdf" , width = 8 , height = 6 )
plot ( my.vars[1:2] , xlab = "PC" , ylab = expression (Q[PC]/F[PC]) , pch = 20 , bty = "n" , xaxt = "n" , yaxt = "n" , ylim = c ( 0 , 60 ) , xlim = c ( 0 , 10 ) , cex = 1.6 )
axis ( 1 , at = seq ( 1,2) )
axis ( 2 , at = seq ( 0 , 60 , 10 ) )
dev.off()



pdf ( "../../../Evolution2014/Figs/HeightPCs1-10.pdf" , width = 8 , height = 6 )
plot ( my.vars , xlab = "PC" , ylab = expression (Q[PC]/F[PC]) , pch = 20 , bty = "n" , xaxt = "n" , yaxt = "n" , ylim = c ( 0 , 60 ) , xlim = c ( 0 , 10 ) , cex = 1.6 )
axis ( 1 , at = seq ( 1,10) )
axis ( 2 , at = seq ( 0 , 60 , 10 ) )
dev.off()


my.pc <- -1*pcas[,2]


plot ( -pcas[,2] , my.heights , pch = 20 , cex = 0.6 , bty = "n" , xlab = "Principal Component 1" , ylab = "Genetic Height Estimates" , ylim = c ( -1 , 1 ) )
abline ( h = 0 , lty = 2 , lwd = 0.7 )
abline ( v = 0 , lty = 2 , lwd = 0.7 )
abline ( a = 0 , b = -sum(my.heights*pcas[,2]) , col = "red" , lwd = 2 )
#abline ( lm ( my.heights ~ my.pc ) )
#abline ( a = 0 , b = -sqrt(qchisq ( 0.95,1)*var) , col = "red")
#abline ( a = 0 , b = 0 , col = "green")




plot ( -pcas[,5] , my.heights , pch = 20 , cex = 0.6 , bty = "n" , xlab = "Principal Component 4" , ylab = "Genetic Height Estimates" , ylim = c ( -1 , 1 ) )
abline ( h = 0 , lty = 2 , lwd = 0.7 )
abline ( v = 0 , lty = 2 , lwd = 0.7 )
abline ( a = 0 , b = -sum(my.heights*pcas[,5]) , col = "red" , lwd = 2 )






pc.eig <- my.pc * eig.vals [ 1, ]
summary(lm ( my.heights ~ my.pc ))
summary(lm ( my.heights ~ pc.eig ))













######## Pigment ##########

setwd("~/Documents/Academics/MultiTraitPolygenicAdaptation/PopData/Pigment")
# temp <- read.delim ( "pigment.genos.all",h=F,sep = " " , stringsAsFactors=F)
# haploid.genos <- temp[,6:ncol(temp)]
# diploid.genos <- matrix ( NA , nrow = nrow ( temp ) , ncol = 2838 )
# #diploid.genos [ , 1:5 ] <- as.matrix(temp [ , 1:5])
# for ( i in 1:2838 ) {
	# diploid.genos [ , i ] <- haploid.genos [ , (i-1)*2+1 ] + haploid.genos [ , (i-1)*2+2 ]
	
# }
# my.data <- cbind ( temp[,1:5],diploid.genos)
# diploid.genos.pigment <- my.data
# new.diploid.genos.pigment<-diploid.genos.pigment[!duplicated ( diploid.genos.pigment$V2),]
# diploid.genos.pigment<-new.diploid.genos.pigment
# diploid.genos.pigment <- diploid.genos.pigment [ order ( diploid.genos.pigment$V2 ) ,  ]
# save ( diploid.genos.pigment , file = "pigment.diploid.genos")
load ( "pigment.diploid.genos" )
pos.chrom <- read.delim ( "pigment.SNP.POS.txt" , h = F )
pos.chrom <- pos.chrom [ pos.chrom$V1 %in% diploid.genos.pigment$V2 , ]
diploid.genos.pigment <- diploid.genos.pigment [ match ( pos.chrom$V1 , diploid.genos.pigment$V2 ) , ]
distances <- matrix ( NA , nrow = nrow ( pos.chrom ) , ncol = nrow ( pos.chrom ) )
for ( i in 1 : nrow ( pos.chrom ) ) {
	for ( j in 1 : nrow ( pos.chrom ) ) {
		if ( pos.chrom$V2 [ i ] == pos.chrom$V2 [ j ] ) {
			distances [ i , j ] <- abs ( pos.chrom$V3 [ i ] - pos.chrom$V3 [ j ] )
		} else  {
			distances [ i , j ] <- -1
		}
	}
}
diploid.genos.pigment <- diploid.genos.pigment [ , c ( 1:5 , indexes.to.keep ) ]
mean.geno <- rowMeans ( diploid.genos.pigment [ , 6 : ncol ( diploid.genos.pigment ) ] )
names ( mean.geno ) <- diploid.genos.pigment$V2
mean.centered <- as.matrix ( diploid.genos.pigment [ , 6 : ncol ( diploid.genos.pigment ) ] - mean.geno )
rownames ( mean.centered ) <- names ( mean.geno )

SNP.var.mat <- cov ( t ( mean.centered ) )
SNP.var.mat [ which ( distances > 200000 | distances < 0 ) ] <- 0
diag ( SNP.var.mat ) <- 2*mean.geno/2 * ( 1 - mean.geno/2 )


# load pigment gwas data
pigment.gwas <- read.delim ( "../../GWAS_Data/pigment.traits.txt" )
my.pigment.gwas <- pigment.gwas [ pigment.gwas$SNP %in% diploid.genos.pigment$V2 , ]




check.vec <- character ( nrow ( my.pigment.gwas ) )
for ( i in 1 : nrow ( my.pigment.gwas ) ) {
	
	if ( (my.pigment.gwas$A1 [ i ] == diploid.genos.pigment$V4 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ] ) & (my.pigment.gwas$A2 [ i ] == diploid.genos.pigment$V5 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ]) ){
		check.vec [ i ] <- T
	} else if ( (my.pigment.gwas$A1 [ i ] == diploid.genos.pigment$V5 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ]) & (my.pigment.gwas$A2 [ i ] == diploid.genos.pigment$V4 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ]) ) {
		check.vec [ i ] <- F
	} else {
		check.vec [ i ] <- "STRAND!!!!!"
	}
}


my.pigment.gwas <- StrandFlip ( my.pigment.gwas , check.vec )



for ( i in 1 : nrow ( my.pigment.gwas ) ) {
	
	if ( (my.pigment.gwas$A1 [ i ] == diploid.genos.pigment$V4 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ]) & (my.pigment.gwas$A2 [ i ] == diploid.genos.pigment$V5 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ]) ){
		next
	} else if ( (my.pigment.gwas$A1 [ i ] == diploid.genos.pigment$V5 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ]) & (my.pigment.gwas$A2 [ i ] == diploid.genos.pigment$V4 [ diploid.genos.pigment$V2 %in% my.pigment.gwas$SNP [ i ] ]) ) {
		temp1 <- my.pigment.gwas$A1 [ i ]
		temp2 <- my.pigment.gwas$A2 [ i ]
		my.pigment.gwas$A1 [ i ] <- temp2
		my.pigment.gwas$A2 [ i ] <- temp1
		my.pigment.gwas$EFF [ i ] <- - my.pigment.gwas$EFF [ i ]
	}
}
save ( my.pigment.gwas , file = "pigment.gwas.als.flipped.Robj")










my.effect.mat <- matrix ( 0 , ncol = length ( unique ( my.pigment.gwas$PHENO ) ) , nrow = length ( unique ( my.pigment.gwas$SNP ) ) )
split.pigment.gwas <- split ( my.pigment.gwas , my.pigment.gwas$PHENO )

for ( i in 1 : length ( split.pigment.gwas ) ) {
	
	my.effect.mat [ diploid.genos.pigment$V2 %in% split.pigment.gwas [[ i ]]$SNP , i ] <- split.pigment.gwas [[ i ]]$EFF
	
}


G.mat <- t ( my.effect.mat ) %*% SNP.var.mat %*% my.effect.mat 
G.sqrt <- eigen ( G.mat )$vec %*% diag ( sqrt ( eigen ( G.mat )$val ) ) %*% t ( eigen ( G.mat )$vec )
phenos <- t ( my.effect.mat ) %*% as.matrix ( diploid.genos.pigment [ , 6:ncol ( diploid.genos.pigment ) ] )
mean.phenos <- rowMeans ( phenos )
mean.centered.phenos <- phenos - mean.phenos


my.var.mats <- list ()
for ( i in 1:10 ) {
	my.var.mats [[ i ]] <- solve ( G.sqrt ) %*% mean.centered.phenos %*% ( pcas[,i+1] / sqrt ( eig.vals [ i , ] ) ) %*% t ( pcas[,i+1] / sqrt ( eig.vals [ i , ] ) ) %*% t ( mean.centered.phenos ) %*% solve ( G.sqrt )
}

my.var.pig <- list ()
for ( j in 1 : 12 ) {
	my.var.pig [[ j ]] <- numeric ( 10 )
	for ( i in 1:10 ) {
		my.var.pig [[ j ]] [ i ] <- solve ( G.mat [ j , j ] ) * mean.centered.phenos [ j , ] %*% ( pcas[,i+1] / sqrt ( eig.vals [ i , ] ) ) * t ( pcas[,i+1] / sqrt ( eig.vals [ i , ] ) ) %*% mean.centered.phenos [ j , ] * solve ( G.mat [ j , j ] )
	}
}


my.eigen.vals <- numeric ( 10 )
for ( i in 1 : 10 ) {
	
	my.eigen.vals [ i ] <- eigen ( my.var.mats [[ i ]] )$va [ 1 ]
	
	
}


library ( MCMCpack )
lapply ( my.var.mats , function ( x ) log ( dwish ( x , length ( indexes.to.keep) , diag ( 12 ) ) ) )




StrandFlip <- function ( gwas , flip.vec ) {
	#recover()
	flip.mat <- matrix ( c ( "A" , "T" , "C" , "G" , "T" , "A" , "G" , "C" ) , ncol = 2 )
	
	for ( i in 1 : length ( flip.vec ) ) {
		if ( flip.vec [ i ] == "STRAND!!!!!" ) {
			gwas [ i , 3 ] <- flip.mat [ flip.mat [ , 1 ] %in% gwas [ i , 3 ] , 2 ]
			gwas [ i , 4 ] <- flip.mat [ flip.mat [ , 1 ] %in% gwas [ i , 4 ] , 2 ]
		}
	}
	
	return ( gwas )
	
}

















######## Histogram cartoon ###########

pop1 <- rnorm ( 1000000 , 1 , 0.4 )
pop2 <- rnorm ( 1000000 , -1 , 0.4 )

setwd ( "/Users/JeremyBerg/Documents/Academics/Evolution2014/Figs" )
pdf ( "two_trait_histo.pdf" , height = 5 , width = 5 )
pdf ( "two_trait_histo_w_lines.pdf" , height = 5 , width = 5 )
pdf ( "two_trait_histo_w_expect.pdf" , height = 5 , width = 5 )
par ( mar =  c ( 1 , 0 , 2 , 1))
plot ( NA , xlim = c ( -2.3 , 2.3 ) , ylim = c ( 0 , 1.1 ) , xaxt = "n" , yaxt = "n" , bty = "n" , xlab = "" , ylab = "")
mtext ( "Frequency" , 2 , -1 )
mtext ( "Phenotype" , 1 , -0.5 )
hist ( pop1 , col = rgb ( 1 , 0 , 0 , 0.5 ) , breaks = 50 , freq = F , add = T , border = rgb ( 1 , 0 , 0 , 0.01 ) )
hist ( pop2 , col = rgb ( 0 , 0 , 1 , 0.5 ) , breaks = 50 , freq = F , add = T , border = rgb ( 0 , 0 , 1 , 0.01 ))
abline ( v = 1 , lty = 2 )
abline ( v = -1 , lty = 2 )

arrows ( -1 , 1.05 , 1 , 1.05 , col = rgb ( 0.3 , 0.6 , 0.3 , 1 ) , angle = 45 , code = 3 , lwd = 5 )
text ( x = 0 , y = 1.1 , expression(V[B]) , col = rgb ( 0.3 , 0.6 , 0.3 , 1 ) , cex = 1.5)

arrows ( -1.4 , 0.4 , -0.6 , 0.4 , col = rgb ( 0 , 0 , 1 , 1 ) , angle = 45 , code = 3 , lwd = 5 )
text ( x = -1.2 , y = 0.55 , expression(V[W]) , col = rgb ( 0 , 0 , 1 , 1 ) , cex = 1.5 )

arrows ( 1.4 , 0.4 , 0.6 , 0.4 , col = rgb ( 1 , 0 , 0 , 1 ) , angle = 45 , code = 3 , lwd = 5 )
text ( x = 0.8 , y = 0.55 , expression(V[W]) , col = rgb ( 1 , 0 , 0 , 1 ) , cex = 1.5 )

arrows ( -0.7 , 0.95 , 0.7 , 0.95 , col = rgb ( 0 , 0 , 0 , 1 ) , angle = 45 , code = 3 , lwd = 5 )
arrows ( -1.7 , 0.3 , -0.3 , 0.3 , col = rgb ( 0 , 0 , 0 , 1 ) , angle = 45 , code = 3 , lwd = 5 )
arrows ( 1.7 , 0.3 , 0.3 , 0.3 , col = rgb ( 0 , 0 , 0 , 1 ) , angle = 45 , code = 3 , lwd = 5 )

dev.off()






par ( mar = c ( 1,1,1,1))
plot ( NA , bty = "n" , xaxt = "n" , yaxt = "n" , ylim = c ( 0 , 1.03 ) , xlim = c ( 0 , 1 ) , ylab = "" , xlab = "")
#lines ( c ( 0.5 , 0.5) , c ( 0.05 , 1) , lwd = 3 )
lines ( c ( 0.2 , 0.5) , c ( 0.05 , 1) , lwd = 3 )
lines ( c ( 0.8 , 0.5) , c ( 0.05 , 1) , lwd = 3 )
text ( x = 0.1 , y = 0.1 , "Present" , cex = 1.5 )
text ( x = 0.1 , y = 1 , "Past" , cex = 1.5 )
text ( x = 0.5 , 1.03 , "Ancestor" , cex = 1.5 )
text ( x = 0.2 , 0 , "Population 1" , cex = 1.5 )
text ( x = 0.8 , 0 , "Population 2" , cex = 1.5 )

norm.1 <- rnorm ( 100 , sd = 0.05 )
norm.2 <- rnorm ( 100 , sd = 0.05 )
norm.3 <- rnorm ( 100 , sd = 0.05 )
norm.4 <- rnorm ( 100 , sd = 0.05 )
# positive
pdf ( "2PopPCPlotPos.pdf" , height = 6 , width = 8 )
pdf ( "2PopPCPlotPosWVb.pdf" , height = 6 , width = 8 )
pdf ( "2PopPCPlotPosWVbFst.pdf" , height = 6 , width = 8 )
par ( mar = c ( 2 , 2 , 1 ,1 ) )
plot ( NA , bty = "n" , xaxt = "n" , yaxt = "n" , ylim = c ( -1 , 1 ) , xlim = c ( -1 , 1 ) , ylab = "" , xlab = "")
abline ( h = 0 , lty = 2 )
abline ( v = 0 , lty = 2 )

points ( 0.8 + norm.1 , 0.2 + norm.2 , pch = 20 )
points ( -0.8 + norm.3 , -0.2 + norm.4 , pch = 20 )
mtext ( "Principal Component 1" , side = 1 , line = 0 )
mtext ( "Phenotype" , side = 2 , line = 0 )
abline ( a = 0 , b = 0.25 , col = "red" , lwd = 2 )

abline ( h = 0.2 , col = rgb ( 0.3 , 0.6 , 0.3 , 1 ) , lty = 3 )
abline ( h = -0.2 , col = rgb ( 0.3 , 0.6 , 0.3 , 1 ) , lty = 3 )
arrows ( -0.4 , 0.2 , -0.4 , -0.2 , code = 3 , col = rgb ( 0.3 , 0.6 , 0.3 , 1 ) , lwd = 3 )
text ( x = -0.3 , y = 0.1 , expression(V[B]) , col = rgb ( 0.3 , 0.6 , 0.3 , 1 ) , cex = 1.5)

abline ( v = 0.8 , lty = 3 )
abline ( v = -0.8 , lty = 3 )
arrows ( -0.8 , 0.35 , 0.8 , 0.35 , code = 3 , lwd = 3 )
text ( x = 0.2 , y = 0.43 , expression(F[ST]) , cex = 1.5)
dev.off ( )


# negative
pdf ( "2PopPCPlotNeg.pdf" , height = 6 , width = 8 )
par ( mar = c ( 2 , 2 , 1 ,1 ) )
plot ( NA , bty = "n" , xaxt = "n" , yaxt = "n" , ylim = c ( -1 , 1 ) , xlim = c ( -1 , 1 ) , ylab = "" , xlab = "")
abline ( h = 0 , lty = 2 )
abline ( v = 0 , lty = 2 )

points ( 0.8 + rnorm ( 100 , sd = 0.05 ) , -0.2 + rnorm ( 100 , sd = 0.05 ) , pch = 20 )
points ( -0.8 + rnorm ( 100 , sd = 0.05 ) , 0.2 + rnorm ( 100 , sd = 0.05 ) , pch = 20 )
mtext ( "Principal Component 1" , side = 1 , line = 0 )
mtext ( "Phenotype" , side = 2 , line = 0 )
abline ( a = 0 , b = -0.25 , col = "red" , lwd = 2 )
dev.off()





