setwd ( "/Users/JeremyBerg/Documents/Academics/MultiTraitPolygenicAdaptation/FinnishSoftware/driftsel/R")
sapply(list.files(pattern="*.R"), source, .GlobalEnv)
setwd("../data")
sapply(list.files(pattern="*.rda"), load, .GlobalEnv)



# genotype data
setwd("~/Documents/Academics/MultiTraitPolygenicAdaptation/PopData/")
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








