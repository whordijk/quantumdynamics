colors<- c( rgb(189,215,231,maxColorValue=555),
            rgb(189,215,231,maxColorValue=255),
            rgb(107,174,214,maxColorValue=555),
            rgb(107,174,214,maxColorValue=255),
            rgb(186,228,179,maxColorValue=555),
            rgb(186,228,179,maxColorValue=255))

#int_single <- read.table("./single_slit.dat")
#int_double <- read.table("./double_slit.dat")
int <- read.table("./intensity.dat")

pdf(file="intensity.pdf",height=4.0,width=5.0)

#plot(int_single$V2,int_single$V1,type="l",col=colors[3],ann=FALSE,xlim=c(0,400),ylim=c(0,100))
#lines(int_double$V2,int_double$V1,col=colors[2])
plot(int$V2,int$V1,type="l",col=colors[3],ann=FALSE,xlim=c(0,400),ylim=c(0,0.05))
title(xlab="y",
        ylab="Intensity")


dev.off()
