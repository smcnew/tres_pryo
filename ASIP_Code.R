d<-read.delim("ASIP/Meth_Data.txt")
d2<-read.delim("ASIP/Align_to_ASIP.txt")

head(d2)
d2$Location<-seq(1,nrow(d2),1)

cpgnum<-66

library(car)

d$sresp1<-d$stress1-d$base1
d$sresp3<-d$stress3-d$base3
d$dresp1<-d$stress1-d$dex1
d$dresp3<-d$stress3-d$dex3
d$l_sresp1<-log(d$sresp1+6)
d$l_sresp3<-log(d$sresp3)
d$l_dresp1<-log(d$stress1)-log(d$dex1)
d$l_dresp3<-log(d$stress3)-log(d$dex3)

assays<-read.delim("ASIP/ASIP_Assay_Locs.txt")

head(assays)
# Make a barplot of percent methylation across the entire region sampled
	pdf("Output_Methylation_Profile_ASIP.pdf",width=9.5,height=4.6)
		par(fig=c(0,1,0,1))
		d3<-d2
		plot(d3$Location,d3$Meth,bty="n",ylim=c(0,100),pch=21,bg="gray70",xaxt="n",ylab="Average methylation %",
			xlab="CpG site number",type="n",las=2,yaxs="i")
		for(i in 1:nrow(d3)){
			lines(rep(d3$Location[i],2),c(d3$Meth[i]+d3$SD[i],d3$Meth[i]-d3$SD[i]),col="gray50")
		}
		points(d3$Location,d3$Meth,pch=21,bg="gray70")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(v=73.5,lty=2,col="red")

## Correlation along sites
	storagec<-rep(NA,cpgnum-1)
	for(i in 1:(cpgnum-1)){
		storagec[i]<-cor(d[,i+9],d[,i+10],use="pairwise.complete.obs")
	}

	colset<-rep(NA,cpgnum-1)
	for(i in 1:length(storagec)){
		ifelse(storagec[i]>0,colset[i]<-"blue",colset[i]<-"red")
		}

	plot(d3$Location[1:(cpgnum-1)],storagec,bty="n",ylim=c(-1,1),pch=21,bg="slateblue",xaxt="n",ylab="Correlation",
		xlab="CpG site number",type="n",las=2,yaxs="i",main="Correlation between adjacent CpG sites")
	abline(h=0,lty=2,col="gray60")
	lines(d3$Location[1:(cpgnum-1)],storagec)
	points(d3$Location[1:(cpgnum-1)],storagec,pch=21,bg=colset)

		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))

	abline(v=73.5,lty=2,col="red")

	yer<- -.75
	len<- .1
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}

## Color treatment
	par(fig=c(0,1,0,1))

	for(i in 1: cpgnum){
			d3$ColCont[i]<-mean(na.omit(subset(d[,i+8],d$Treat1=="Control")))
			d3$ColDull[i]<-mean(na.omit(subset(d[,i+8],d$Treat1=="Dull")))
			d3$StrCont[i]<-mean(na.omit(subset(d[,i+8],d$Treat2=="Control")))
			d3$Stress[i]<-mean(na.omit(subset(d[,i+8],d$Treat2=="Stress")))
	}

	plot(d3$Location,d3$ColCont,bty="n",ylim=c(0,118),pch=21,bg="white",xaxt="n",ylab="Average methylation %",
			xlab="CpG site number in relation to ATG",type="n",las=2,yaxs="i")
		#for(i in 1:nrow(d3)){
		#	lines(rep(d3$Location[i],2),c(d3$Meth[i]+d3$SD[i],d3$Meth[i]-d3$SD[i]),col="gray50")
		#}
		points(d3$Location,d3$ColCont,pch=21,bg="white")
		points(d3$Location,d3$ColDull,pch=21,bg="gray55")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		lines(rep(-7.5,2),c(0,99),lty=2,col="red")
		abline(v=73.5,lty=2,col="red")

		legend("topleft",c("Color Control","Color Dull"),pch=c(21,21),pt.bg=c("white","gray55"),bty="n")


	plot(d3$Location,d3$StrCont,bty="n",ylim=c(0,118),pch=21,bg="white",xaxt="n",ylab="Average methylation %",
			xlab="CpG site number in relation to ATG",type="n",las=2,yaxs="i")
		#for(i in 1:nrow(d3)){
		#	lines(rep(d3$Location[i],2),c(d3$Meth[i]+d3$SD[i],d3$Meth[i]-d3$SD[i]),col="gray50")
		#}
		points(d3$Location,d3$StrCont,pch=21,bg="slateblue")
		points(d3$Location,d3$Stress,pch=21,bg="orange")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(v=73.5,lty=2,col="red")

		legend("topleft",c("Stress Control","Stress"),pch=c(21,21),pt.bg=c("slateblue","orange"),bty="n")

## Color treatment effect size
	par(fig=c(0,1,0,1))

	for(i in 1: cpgnum){
			d3$ColCont[i]<-mean(na.omit(subset(scale(d[,i+8]),d$Treat1=="Control")))
			d3$ColDull[i]<-mean(na.omit(subset(scale(d[,i+8]),d$Treat1=="Dull")))
			d3$StrCont[i]<-mean(na.omit(subset(scale(d[,i+8]),d$Treat2=="Control")))
			d3$Stress[i]<-mean(na.omit(subset(scale(d[,i+8]),d$Treat2=="Stress")))
	}



	d3$col.diff<-d3$ColCont-d3$ColDull
	d3$str.diff<-d3$StrCont-d3$Stress

	for(i in 1:nrow(d3)){
		ifelse(d3$col.diff[i]>0,d3$color[i]<-"slateblue",d3$color[i]<-"orange")
		}

	plot(d3$Location,d3$col.diff,bty="n",ylim=c(-2,2),pch=21,bg="slateblue",xaxt="n",ylab="Standardized difference",
		xlab="CpG site number",type="n",las=2,yaxs="i",main="Color treatment effects")
	abline(h=0,lty=2,col="gray60")
	lines(d3$Location,d3$col.diff)
	points(d3$Location,d3$col.diff,pch=21,bg=d3$color)

		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))

	text(-1,1,"Control ++",col="slateblue",srt=90)
	text(-1,-1,"Color ++",col="orange",srt=90)

		abline(v=73.5,lty=2,col="red")

	yer<- -1.25
	len<- .1
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}

	for(i in 1:nrow(d3)){
		ifelse(d3$str.diff[i]>0,d3$color[i]<-"slateblue",d3$color[i]<-"orange")
		}

	plot(d3$Location,d3$str.diff,bty="n",ylim=c(-2,2),pch=21,bg="slateblue",xaxt="n",ylab="Standardized difference",
		xlab="CpG site number",type="n",las=2,yaxs="i",main="Stress treatment effects")
	abline(h=0,lty=2,col="gray60")
	lines(d3$Location,d3$str.diff)
	points(d3$Location,d3$str.diff,pch=21,bg=d3$color)

		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))

	text(-1,1,"Control ++",col="slateblue",srt=90)
	text(-1,-1,"Stress ++",col="orange",srt=90)

		yer<- -1.25
	len<- .1
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}


## Make plot of one way t-tests for each predictor

	## Color manipulation treatment
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-t.test(subset(d[,i+8],d$Treat1=="Control"),subset(d[,i+8],d$Treat1=="Dull"))
			storagel[i]<- -log(test$p.value,base=10)
			storagep[i]<-test$p.value
			}
		plot(d3$Location,storagel,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number",bty="n",main="Color Manipulation T-Tests")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

			yer<- 4
	len<- -.2
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}

	## Stress manipulation treatment
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-t.test(subset(d[,i+8],d$Treat2=="Control"),subset(d[,i+8],d$Treat2=="Stress"))
			storagel[i]<- -log(test$p.value,base=10)
			storagep[i]<-test$p.value
			}
		plot(d3$Location,storagel,pch=21,bg="lightblue",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress Manipulation T-Tests")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

				yer<- 4
	len<- -.2
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}

	## Age
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-t.test(subset(d[,i+8],d$Age=="ASY"),subset(d[,i+8],d$Age=="SY"))
			storagel[i]<- -log(test$p.value,base=10)
			storagep[i]<-test$p.value
			}
		plot(d3$Location,storagel,pch=21,bg="yellow",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Age T-Tests")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
		abline(v=73.5,lty=2,col="red")

					yer<- 4
	len<- -.2
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}

	## breast brightness regression
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		storagee<-rep(NA,cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$b.bright)
			test2<-lm(scale(d[,i+8])~scale(d$b.bright))
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
			storagee[i]<-coef(test2)[2]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Breast Brightness Regression")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")


					colset<-rep(NA,length(storagee))
		for(i in 1:length(storagee)){
			ifelse(storagee[i]>0,colset[i]<-"slateblue",colset[i]<-"orange")
		}
		plot(d3$Location,storagee,pch=21,bg=colset,ylim=c(-2,2),ylab="Standardized difference",xlab="CpG number",
			xaxt="n",bty="n",main="Stress 3 Effects",las=2,type="n",yaxs="i")
		lines(d3$Location,storagee)
		points(d3$Location,storagee,pch=21,bg=colset)
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(v=73.5,lty=2,col="red")
		abline(h=0,lty=2,col="gray60")

		text(-2,1,"Higher SC ++",col="slateblue",srt=90)
		text(-2,-1,"Lower SC ++",col="orange",srt=90)

						yer<- -1
	len<- .1
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}

	## base cort 1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$base1)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Base Cort 1")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## stress cort 1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$stress1)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress Cort 1")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## dex cort 1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$dex1)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Dex Cort 1")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## base cort 2
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$base2)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Base Cort 2")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## base cort 3
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$base3)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Base Cort 3")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## stress cort 3
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		storagee<-rep(NA,cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$stress3)
			test2<-lm(scale(d[,i+8])~scale(d$stress3))
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
			storagee[i]<-coef(test2)[2]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress Cort 3")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}

				colset<-rep(NA,length(storagee))
		for(i in 1:length(storagee)){
			ifelse(storagee[i]>0,colset[i]<-"slateblue",colset[i]<-"orange")
		}
		plot(d3$Location,storagee,pch=21,bg=colset,ylim=c(-2,2),ylab="Standardized difference",xlab="CpG number",
			xaxt="n",bty="n",main="Stress 3 Effects",las=2,type="n",yaxs="i")
		lines(d3$Location,storagee)
		points(d3$Location,storagee,pch=21,bg=colset)
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(v=73.5,lty=2,col="red")
		abline(h=0,lty=2,col="gray60")

		text(-2,1,"Higher SC ++",col="slateblue",srt=90)
		text(-2,-1,"Lower SC ++",col="orange",srt=90)

						yer<- -1
	len<- .1
	for(i in 1:nrow(assays)){
		arrows(assays$Start1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		arrows(assays$End1[i],yer,y1=yer+len,angle=30,col="black",length=.05)
		lines(c(assays$Start1[i],assays$End1[i]),rep(yer,2))
	}


	## dex cort 3
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$dex3)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Dex Cort 3")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## sresp1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$sresp1)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress response 1")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## sresp1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$sresp3)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress response 3")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## log sresp1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$l_sresp1)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Log Stress response 1")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## log sresp3
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$l_sresp3)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Log Stress response 3")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## dresp1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$dresp1)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Dex response 1")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## dresp3
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$dresp3)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Dex response 3")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## log dresp1
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$l_dresp1)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5.5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Log Dex response 1")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

	## log dresp3
		storagel<-rep(NA, cpgnum)
		storagep<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-lm(d[,i+8]~d$l_dresp3)
			storagel[i]<- -log(summary(test)$coefficients[8],base=10)
			storagep[i]<-summary(test)$coefficients[8]
				}
		plot(d3$Location,storagel,pch=21,bg="gray50",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Log Dex response 3")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

## Run as more complicated models

	# main effects base, stress, dex capture 1
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		storagel3<-rep(NA, cpgnum)

		for(i in 1: cpgnum){
			test<-summary(lm(d[,i+8]~d$base1+d$stress1+d$dex1))
			storagel1[i]<- -log(test$coefficients[2,4],base=10)
			storagel2[i]<- -log(test$coefficients[3,4],base=10)
			storagel3[i]<- -log(test$coefficients[4,4],base=10)
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4.5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress series 1")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		points(d3$Location-0.1,storagel3,pch=21,bg="violet")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Base","Stress","Dex"),pch=21,
			pt.bg=c("orange","lightblue","violet"),bty="n")

	# Interaction between base, stress, dex
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		storagel3<-rep(NA, cpgnum)
		storagel4<-rep(NA, cpgnum)
		storagel5<-rep(NA, cpgnum)
		storagel6<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-summary(lm(d[,i+8]~d$base1*d$stress1+d$base1*d$dex1+d$stress1*d$dex1))
			storagel1[i]<- -log(test$coefficients[2,4],base=10)
			storagel2[i]<- -log(test$coefficients[3,4],base=10)
			storagel3[i]<- -log(test$coefficients[4,4],base=10)
			storagel4[i]<- -log(test$coefficients[5,4],base=10)
			storagel5[i]<- -log(test$coefficients[6,4],base=10)
			storagel6[i]<- -log(test$coefficients[7,4],base=10)
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress series 1 with interactions")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		points(d3$Location-0.1,storagel3,pch=21,bg="violet")
		points(d3$Location+0.1,storagel4,pch=21,bg="gray")
		points(d3$Location-0.1,storagel5,pch=21,bg="white")
		points(d3$Location-0.1,storagel6,pch=21,bg="pink")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Base","Stress","Dex","Base*Stress","Base*Dex","Stress*Dex"),pch=21,
			pt.bg=c("orange","lightblue","violet","gray","white","pink"),bty="n")

	# main effects base, stress, dex capture 2
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		storagel3<-rep(NA, cpgnum)

		for(i in 1: cpgnum){
			test<-summary(lm(d[,i+8]~d$base3+d$stress3+d$dex3))
			storagel1[i]<- -log(test$coefficients[2,4],base=10)
			storagel2[i]<- -log(test$coefficients[3,4],base=10)
			storagel3[i]<- -log(test$coefficients[4,4],base=10)
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4.5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress series 3")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		points(d3$Location-0.1,storagel3,pch=21,bg="violet")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Base","Stress","Dex"),pch=21,
			pt.bg=c("orange","lightblue","violet"),bty="n")

	# Interaction between base, stress, dex capture 2
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		storagel3<-rep(NA, cpgnum)
		storagel4<-rep(NA, cpgnum)
		storagel5<-rep(NA, cpgnum)
		storagel6<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-summary(lm(d[,i+8]~d$base3*d$stress3+d$base3*d$dex3+d$stress3*d$dex3))
			storagel1[i]<- -log(test$coefficients[2,4],base=10)
			storagel2[i]<- -log(test$coefficients[3,4],base=10)
			storagel3[i]<- -log(test$coefficients[4,4],base=10)
			storagel4[i]<- -log(test$coefficients[5,4],base=10)
			storagel5[i]<- -log(test$coefficients[6,4],base=10)
			storagel6[i]<- -log(test$coefficients[7,4],base=10)
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4.5),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Stress series 3 with interactions")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		points(d3$Location-0.1,storagel3,pch=21,bg="violet")
		points(d3$Location+0.1,storagel4,pch=21,bg="gray")
		points(d3$Location-0.1,storagel5,pch=21,bg="white")
		points(d3$Location-0.1,storagel6,pch=21,bg="pink")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Base","Stress","Dex","Base*Stress","Base*Dex","Stress*Dex"),pch=21,
			pt.bg=c("orange","lightblue","violet","gray","white","pink"),bty="n")

	# Two way ANOVA with the two treatments
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		two.way.mods<-list()
		two.way.mods2<-list()
		for(i in 1: cpgnum){
			test<-Anova(lm(d[,i+8]~d$Treat1+d$Treat2),type="III")
			test2<-summary(lm(d[,i+8]~d$Treat1+d$Treat2))
			storagel1[i]<- -log(test["Pr(>F)"][2,1],base=10)
			storagel2[i]<- -log(test["Pr(>F)"][3,1],base=10)
			two.way.mods[[i]]<-test
			two.way.mods2[[i]]<-test2
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Two-way ANOVA Color + Stress Manipulations")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Color treatment","Stress treatment"),pch=21,pt.bg=c("orange","lightblue"),bty="n")


	## Interaction between the two treatments
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		storagel3<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-Anova(lm(d[,i+8]~d$Treat1*d$Treat2),type="III")
			storagel1[i]<- -log(test["Pr(>F)"][2,1],base=10)
			storagel2[i]<- -log(test["Pr(>F)"][3,1],base=10)
			storagel3[i]<- -log(test["Pr(>F)"][4,1],base=10)
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Two-way ANOVA Color*Stress Manipulations")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		points(d3$Location-0.1,storagel3,pch=21,bg="violet")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Color treatment","Stress treatment","Color*Stress"),pch=21,
			pt.bg=c("orange","lightblue","violet"),bty="n")

	## Interaction between the two treatments
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		storagel3<-rep(NA, cpgnum)
		storagel4<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-Anova(lm(d[,i+8]~d$Treat1*d$Treat2+d$b.bright),type="III")
			storagel1[i]<- -log(test["Pr(>F)"][2,1],base=10)
			storagel2[i]<- -log(test["Pr(>F)"][3,1],base=10)
			storagel3[i]<- -log(test["Pr(>F)"][5,1],base=10)
			storagel4[i]<- -log(test["Pr(>F)"][4,1],base=10)
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Two-way ANOVA Color*Stress + Breast Brightness")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		points(d3$Location-0.1,storagel3,pch=21,bg="violet")
		points(d3$Location,storagel4,pch=21,bg="gray55")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Color treatment","Stress treatment","Color*Stress","Breast brightness"),pch=21,
			pt.bg=c("orange","lightblue","violet","gray55"),bty="n")

	## Interaction between the two treatments
		storagel1<-rep(NA, cpgnum)
		storagel2<-rep(NA, cpgnum)
		storagel3<-rep(NA, cpgnum)
		storagel4<-rep(NA, cpgnum)
		storagel5<-rep(NA, cpgnum)
		storagel6<-rep(NA, cpgnum)
		for(i in 1: cpgnum){
			test<-Anova(lm(d[,i+8]~d$Treat1*d$Treat2+d$b.bright*d$Treat1+d$b.bright*d$Treat2),type="III")
			storagel1[i]<- -log(test["Pr(>F)"][2,1],base=10)
			storagel2[i]<- -log(test["Pr(>F)"][3,1],base=10)
			storagel3[i]<- -log(test["Pr(>F)"][5,1],base=10)
			storagel4[i]<- -log(test["Pr(>F)"][4,1],base=10)
			storagel5[i]<- -log(test["Pr(>F)"][6,1],base=10)
			storagel6[i]<- -log(test["Pr(>F)"][7,1],base=10)
				}
		plot(d3$Location,storagel1,pch=21,bg="orange",ylab=expression('-log'[10]*'(P)'),las=2,ylim=c(0,4),xaxt="n",
			xlab="CpG site number in relation to ATG",bty="n",main="Two-way ANOVA Color*Stress + Breast*Treatments")
		points(d3$Location+0.1,storagel2,pch=21,bg="lightblue")
		points(d3$Location-0.1,storagel3,pch=21,bg="violet")
		points(d3$Location,storagel4,pch=21,bg="gray55")
		points(d3$Location+.1,storagel5,pch=21,bg="lightgreen")
		points(d3$Location-.1,storagel6,pch=21,bg="white")
		axis(1,at=c(-300,-250,-200,-150,-100,-50,0,50,100,150,200,250))
		abline(h=1.30103,lty=2,col="gray55")
			text(-276,1.2,"P = 0.05",pos=3,col="gray35",cex=0.8)
		abline(h=2,lty=2,col="gray55")
			text(-276,1.9,"P = 0.01",pos=3,col="gray35",cex=0.8)
		abline(h=3,lty=2,col="gray55")
			text(-276,2.9,"P = 0.001",pos=3,col="gray35",cex=0.8)
			abline(v=73.5,lty=2,col="red")

		legend("topright",c("Color treatment","Stress treatment","Color*Stress","Breast brightness",
			"Breast*Color","Breast*Stress"),pch=21,
			pt.bg=c("orange","lightblue","violet","gray55","lightgreen","white"),bty="n")

	dev.off()

	## random

		storagep<-rep(NA,104)
		storagep2<-rep(NA,104)
		storagep3<-rep(NA,104)
		storagep4<-rep(NA,104)
		rr<-rep(NA,1000)
		rr2<-rep(NA,1000)
		rr3<-rep(NA,1000)
		rr4<-rep(NA,1000)
		for(k in 1:1000){
				preds<-c(rep("A",18),rep("B",18))
				preds2<-c(rep("A",18),rep("B",18))
				set.seed(109)
				rando<-runif(nrow(d),0,1)
				ord<-as.data.frame(cbind(preds,rando))
				ord<-ord[order(ord$rando),]
				d$rando<-ord$preds
				rando<-runif(nrow(d),0,1)
				ord<-as.data.frame(cbind(preds2,rando))
				ord<-ord[order(ord$rando),]
				d$rando2<-ord$preds2
			for(i in 1:104){
				test<-t.test(subset(d[,i+8],d$rando=="A"),subset(d[,i+8],d$rando=="B"))
				test2<-Anova(lm(d[,i+8]~d$rando*d$rando2),type="III")
				storagep[i]<-test$p.value
				storagep2[i]<-test2["Pr(>F)"][2,1]
				storagep3[i]<-test2["Pr(>F)"][3,1]
				storagep4[i]<-test2["Pr(>F)"][4,1]
			}
				rr[k]<-length(subset(storagep,storagep<0.01))
				rr2[k]<-length(subset(storagep2,storagep2<0.01))
				rr3[k]<-length(subset(storagep3,storagep3<0.01))
				rr4[k]<-length(subset(storagep4,storagep4<0.01))
				print(k)
		}


pdf("Output_rando_treat.pdf",width=9,height=4)
		par(mfrow=c(1,3))
		hist(rr,breaks=seq(0,22,1),yaxs="i",yaxt="n",ylim=c(0,1000),xaxs="i",xlab="Number of CpG with P < 0.01",
			ylab="Percentage of permutations",main="One-Way ANOVA")
		axis(2,at=seq(0,1000,100),labels=c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"),las=2)
		points(0,100,pch=21,bg="orange")
		points(2,150,pc=21,bg="lightblue")
		points(0,100,pch=21,bg="gray55")
		abline(v=quantile(rr,0.95),lty=2,col="red")


		hist(c(rr2,rr3),breaks=seq(0,25,1),yaxs="i",yaxt="n",ylim=c(0,2000),xaxs="i",xlab="Number of CpG with P < 0.01",
			ylab="Percentage of permutations",main="Two-Way ANOVA Main Effects")
		axis(2,at=seq(0,2000,200),labels=c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"),las=2)
		points(2,200,pch=21,bg="black")
		abline(v=quantile(c(rr2,rr3),0.95),lty=2,col="red")


		hist(rr4,breaks=seq(0,20,1),yaxs="i",yaxt="n",ylim=c(0,1000),xaxs="i",xlab="Number of CpG with P < 0.01",
			ylab="Percentage of permutations",main="Two-Way ANOVA Interaction")
		axis(2,at=seq(0,1000,100),labels=c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"),las=2)
		points(1,100,pch=21,bg="violet")
		abline(v=quantile(rr4,.95),lty=2,col="red")
	dev.off()

## random cort

	main.p<-c(storagel1,storagel2,storagel3)
	inter.p<-c(storagel4,storagel5,storagel6)
	cycles<-1000
	main.p05<-rep(NA,cycles)
	main.p01<-rep(NA,cycles)
	main.p001<-rep(NA,cycles)
	inter.p05<-rep(NA,cycles)
	inter.p01<-rep(NA,cycles)
	inter.p001<-rep(NA,cycles)
	storagel1.r<-rep(NA,104)
	storagel2.r<-rep(NA,104)
	storagel3.r<-rep(NA,104)
	storagel4.r<-rep(NA,104)
	storagel5.r<-rep(NA,104)
	storagel6.r<-rep(NA,104)

	set.seed(99)
	for(i in 1:cycles){
		d.rand<-d[,c("base3","stress3","dex3")]
		d.rand$sort<-runif(nrow(d.rand))
		d.rand<-d.rand[order(d.rand$sort),]
		colnames(d.rand)<-c("base3.r","stress3.r","dex3.r","sort")
		d2<-d
		d2<-cbind(d2,d.rand)
		for(k in 1:104){
			test<-summary(lm(d2[,k+8]~d2$base3.r*d2$stress3.r+d2$base3.r*d2$dex3.r+d2$stress3.r*d2$dex3.r))
					storagel1.r[k]<- -log(test$coefficients[2,4],base=10)
					storagel2.r[k]<- -log(test$coefficients[3,4],base=10)
					storagel3.r[k]<- -log(test$coefficients[4,4],base=10)
					storagel4.r[k]<- -log(test$coefficients[5,4],base=10)
					storagel5.r[k]<- -log(test$coefficients[6,4],base=10)
					storagel6.r[k]<- -log(test$coefficients[7,4],base=10)
				}
		main.r<-c(storagel1.r,storagel2.r,storagel3.r)
		inter.r<-c(storagel4.r,storagel5.r,storagel6.r)
		main.p05[i]<-length(subset(main.r,main.r> -log(0.05,base=10)))
		main.p01[i]<-length(subset(main.r,main.r> -log(0.01,base=10)))
		main.p001[i]<-length(subset(main.r,main.r> -log(0.001,base=10)))
		inter.p05[i]<-length(subset(inter.r,inter.r> -log(0.05,base=10)))
		inter.p01[i]<-length(subset(inter.r,inter.r> -log(0.01,base=10)))
		inter.p001[i]<-length(subset(inter.r,inter.r> -log(0.001,base=10)))
		print(paste(i," of ",cycles,sep=""))
		}


		main.obs.05<-length(subset(main.p,main.p<0.05))
		main.obs.01<-length(subset(main.p,main.p<0.01))
		main.obs.001<-length(subset(main.p,main.p<0.001))
		inter.obs.05<-length(subset(inter.p,inter.p<0.05))
		inter.obs.01<-length(subset(inter.p,inter.p<0.01))
		inter.obs.001<-length(subset(inter.p,inter.p<0.001))


		text.y<-650

pdf("Output_rando_cort.pdf",height=8,width=9)
		par(mfrow=c(2,3))
		hist(main.p05,xlab="P<0.05",ylab="",main="Main effects",ylim=c(0,700))
		abline(v=main.obs.05,lty=2,col="red")
		p1<-round((length(subset(main.p05,main.p05>main.obs.05))+.000001)/cycles,3)
		text(main.obs.05,text.y,paste("P = ",p1,sep=""))

		hist(main.p01,xlab="P<0.01",ylab="",main="Main effects",ylim=c(0,700))
		abline(v=main.obs.01,lty=2,col="red")
		p2<-round((length(subset(main.p01,main.p01>main.obs.01))+.000001)/cycles,3)
		text(main.obs.01,text.y,paste("P = ",p2,sep=""))

		hist(main.p001,xlab="P<0.001",ylab="",main="Main effects",ylim=c(0,700))
		abline(v=main.obs.001,lty=2,col="red")
		p3<-round((length(subset(main.p001,main.p001>main.obs.001))+.000001)/cycles,3)
		text(main.obs.001,text.y,paste("P = ",p3,sep=""))

		hist(inter.p05,xlab="P<0.05",ylab="",main="Interaction effects",ylim=c(0,700))
		abline(v=inter.obs.05,lty=2,col="red")
		p4<-round((length(subset(inter.p05,inter.p05>inter.obs.05))+.000001)/cycles,3)
		text(inter.obs.05,text.y,paste("P = ",p4,sep=""))

		hist(inter.p01,xlab="P<0.01",ylab="",main="Interaction effects",ylim=c(0,700))
		abline(v=inter.obs.01,lty=2,col="red")
		p5<-round((length(subset(inter.p01,inter.p01>inter.obs.01))+.000001)/cycles,3)
		text(inter.obs.01,text.y,paste("P = ",p5,sep=""))

		hist(inter.p001,xlab="P<0.001",ylab="",main="Interaction effects",ylim=c(0,700))
		abline(v=inter.obs.001,lty=2,col="red")
		p6<-round((length(subset(inter.p001,inter.p001>inter.obs.001))+.000001)/cycles,3)
		text(inter.obs.001,text.y,paste("P = ",p6,sep=""))
dev.off()

