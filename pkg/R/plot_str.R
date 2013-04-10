plot_str <-
function(loc=NULL,list_jobs,x_width,y_height){
	list_col=c('darkviolet','magenta','olivedrab3','yellow','purple2','sienna1','wheat1','orangered1','springgreen','navy','seagreen3')
#pdf(file=name_pdf,width=x_width,height=y_height)
	
	for (j in 1:length(list_jobs)){
		job=list_jobs[j]
		datname=paste(loc,'results_job_',job,'_q',sep='')
		dat=as.matrix(read.table(datname))
		npop=length(dat[1,])-2
		probs=NULL
		name_pdf=paste(loc,'job_',job,'.pdf',sep='')
		pdf(file=name_pdf,width=x_width,height=y_height)
		
		
		for (k in 3:(npop+2)){
			#p=as.numeric(dat[,k])
			probs=cbind(probs,as.numeric(dat[,k]))
		}
		
############# get population limits ###########
		pops=dat[,2]

		lpop=unique(pops)
		lim=c(0)
		nr=1:length(pops)
		for (i in 1:length(lpop)){
			fff=lpop[i]
			sub=subset(nr,pops==fff)
			l=max(sub)
			lim=c(lim,l)
			
		}
		
###############################################
		if (npop<=length(list_col)){
			barplot(t(probs),border = NA,col=list_col[1:npop],space=0,axis.lty=1)
		}
		
		if (npop>length(list_col)){
			n=npop-length(list_col)
			extra_col=seq(from=10,to=657,length.out=n)
			ecol=c(list_col,colors()[extra_col])
			barplot(t(probs),border = NA,col=ecol,space=0)
		}
		
		abline(v=lim)
		for (l in 2:length(lim)){
			l1=lim[l-1]
			l2=lim[l]
#print(c(l1,l2,mean(c(l1,l2))))
			mtext(lpop[(l-1)],1,at=mean(c(l1,l2)))
		}
		dev.off()
		
			
		
	}




}
