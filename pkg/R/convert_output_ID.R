convert_IDs <-
function(ID1=NULL,out_nam=NULL,rev_convert=NULL,printqhat=NULL,k=NULL){
		
	if (printqhat==1){
		
		### convert _q file ###
		res_name=paste(out_nam,'_q',sep='')  # read _q output file
		res=as.matrix(read.table(res_name))
		res[,2]=ID1							# change column of IDs by the old one
		write(t(res),ncolumns=length(res[1,]),file=res_name)  # write the _q file back in the result folder .... simple
		##########
		
	}
	
		

	###### convert _f file .... more complicated
	n_max=k
	ncol=(10*n_max)+8  # maximum number of column in the _f file (Otherwise R would arbitrarily fix column number and split large lines into two)
	c_name=numeric(ncol)
	for (i in 1:ncol){
		c_name[i]=paste('V',as.character(i),sep='')  # create dummy column names to force R to read a matrix with (10*n_max)+8 columns
	}												# this avoids R to split one line into two arbitrarily when fixes the number of column of the matrix
	
	res_name=paste(out_nam,'_f',sep='')
	res_scan=scan(res_name,sep='\n',what='character',blank.lines.skip = F)
	res_tab=as.matrix(read.table(res_name,fill=T,sep='',strip.white=F,col.names=c_name))
	n_pop=length(rev_convert)
	ans1='Given'
	ans2='(%Miss)'
	
	## read along _f file 
	
	write(res_scan[1],file=res_name)  # write 1st line of the result file
	i=1
	done='F'
	while (done=='F'){
		i=i+1
		l=strsplit(res_scan[i],' ')[[1]]
		if (!(ans1 %in% l) & !(ans2 %in% l)) write(res_scan[i],file=res_name,append=T)  # if generic line -> write it on the _f file
		if ((ans1 %in% l)){  # if ans1 in the current line: we found the proportion of membership table 
			write(res_scan[i],file=res_name,append=T)
			i=i+1
			l=strsplit(res_scan[i],' ')[[1]]
			write(res_scan[i],file=res_name,append=T)
			###### now change the IDs in the membership table 
			i=i+1
			for (j in (i+1):(i+n_pop)){
				l1=strsplit(res_scan[j],':')[[1]]
				l2=as.character(as.numeric(l1[1]))
				l3=rev_convert[[l2]]
				l4=paste('  ',l3,':',l1[2],sep='')
				write(l4,file=res_name,append=T)
			}
			i=i+n_pop
				
			
			
		}
		l=strsplit(res_scan[i],' ')[[1]]
		if (ans2 %in% l){  # if ans2 in the current line: we found the individual assignment table 
		
			write(res_scan[i],file=res_name,append=T)	
			
			###### now look for the assignment table in the res_tab
			for (j in 1:length(res_tab[,1])){
				if (ans2 %in% res_tab[j,]){  # find the assignment table in res_tab
					
					ref_line=j
					sdat=res_tab[(ref_line+1):(ref_line+length(ID1)),]  # extract teh assignmnet table in "sdat" object
					sdat[,4]=ID1 # change the ID column by the old IDs 
					write(t(sdat),ncolumn=length(sdat[1,]),file=res_name,append=T)
						
				}
			}
			i=i+length(ID1)+1	
		}
		if (i==length(res_scan)) done='T'
		
		
	}
	#########
}



