
####### Copyright Francois besnier 24/01/2013  (The program is distributed under the terms of the GNU General Public License)
####### francois.besnier@imr.no

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


parallel_structure <-
function(joblist=NULL,n_cpu=NULL,structure_path=Mac_path,infile=NULL,outpath=NULL,numinds=NULL,numloci=NULL,
		plot_output=1,label=1,popdata=1,popflag=0,locdata=0,phenotypes=0,markernames=0,mapdist=0,onerowperind=0,phaseinfo=0,
		recessivealleles=0,phased=0,extracol=0,missing=-9,ploidy=2,noadmix=0,linkage=0,usepopinfo=0,locprior=0,
		inferalpha=1,alpha=1.0,popalphas=0,unifprioralpha=1,alphamax=10.0,alphapropsd=0.025,freqscorr=1,onefst=0,
		fpriormean=0.01,fpriorsd=0.05,inferlambda=0,lambda=1.0,computeprob=1,pfromflagonly=0,ancestdist=0,
		startatpopinfo=0,metrofreq=10,updatefreq=1,printqhat=0){
			
	
		require(parallel)
			
	# create a list of global parameters common to all jobs. To be shared with slave nodes
	GlobPar=list(numloci=numloci,label=label,popdata=popdata,popflag=popflag,locdata=locdata,phenotypes=phenotypes,markernames=markernames,mapdist=mapdist,onerowperind=onerowperind,phaseinfo=phaseinfo,recessivealleles=recessivealleles,phased=phased,extracol=extracol,missing=missing,ploidy=ploidy,noadmix=noadmix,linkage=linkage,usepopinfo=usepopinfo,locprior=locprior,inferalpha=inferalpha,alpha=alpha,popalphas=popalphas,unifprioralpha=unifprioralpha,alphamax=alphamax,alphapropsd=alphapropsd,freqscorr=freqscorr,onefst=onefst,fpriormean=fpriormean,fpriorsd=fpriorsd,inferlambda=inferlambda,lambda=lambda,computeprob=computeprob,pfromflagonly=pfromflagonly,ancestdist=ancestdist,startatpopinfo=startatpopinfo,metrofreq=metrofreq,updatefreq=updatefreq,printqhat=printqhat)

	mes=paste('starting work at ',Sys.time(),sep='')
	print(mes)
			
	# always write a file called 'extraparams'  actually empty as all paramaters are in dynamicly generated parameter files for each job
	# structure crashes if no "extramarams" file present in working directory...
	write('',file='extraparams')   
			
			
	
	###### Check data input ######
	###############################
			
	data_name=infile
	n_loci=numloci
	
	###############################


			if (markernames==1){
				mark_line=as.matrix(read.table(data_name,nrows=1))
				dat=as.matrix(read.table(data_name,skip=1))
			}
			
			if (markernames==0){
				mark_line=NULL
				dat=as.matrix(read.table(data_name))
			}
# list_all_pop=unique(dat[,1])  # 08042013
			if (label==1) list_all_pop=unique(dat[,2])
			if (label==0) list_all_pop=unique(dat[,1])
			Npop=length(list_all_pop)

						
			if (onerowperind==0) nind=(length(dat[,1]))/2
			if (onerowperind==1) nind=length(dat[,1])
			
			
	if (nind!=numinds) {
		m=paste('Error: expected ',numinds,' individuals and found ', nind,' please check data file')
		stop(m)
	}
				
	###############################
						
	#### Import job list ##########
						
	job_list=as.matrix(read.table(joblist))
	n_tasks=length(job_list[,1])
	tasks=vector('list')
	for (i in 1:n_tasks){
		job=job_list[i,]
		tasks[[i]]=list(id=as.character(job[1]),job_pop=job[2],k=as.character(job[3]),burnin=as.character(job[4]),iter=as.character(job[5]),path=structure_path)
		
	}
			
	###### if full pairwise matrix in joblist ##########
	for (i in 1:length(tasks)){
		if (tasks[[i]]$job_pop=='pairwise.matrix'){  # if pairwise matrix 
														##  -1 remove the task from the list
														## -2 makes extra tasks for each pair of population in the datafile
			Mat_id=tasks[[i]]$id
			ki=tasks[[i]]$k
			burnini=tasks[[i]]$burnin
			iteri=tasks[[i]]$iter
			tasks[[i]]=NULL
			for (j in 1:(length(list_all_pop)-1)){
				for (k in (j+1):length(list_all_pop)){
					popi1=paste(list_all_pop[j],',',list_all_pop[k],sep='')
					popi2=paste(list_all_pop[j],'_',list_all_pop[k],sep='')
					l=length(tasks)
					jobid=paste(Mat_id,'_',popi2,sep='')
					tasks[[l+1]]=list(id=jobid,job_pop=popi1,k=ki,burnin=burnini,iter=iteri,path=structure_path)
				}
			}
		}
	}
	ALL_tasks=tasks
			
			
			
	slave_fun_unix=function(job){
#Mac_path="/Applications/Structure.app/Contents/Resources/Java/bin/"	
		structure_path=job$path
		id=job$id
		job_pop=job$job_pop   # list populaiton for this job
		k=as.character(job$k)
		burnin=as.character(job$burnin)
		iter=as.character(job$iter)
		
		########################
		temp_nam=paste('temp_list_job_',as.character(id),sep='')
		write(job_pop,file=temp_nam)
		job_pop=as.matrix(read.table(temp_nam,sep=','))
		instr=paste('rm ',temp_nam,sep='')
		system(instr)
		
		#### create a data subset for the poulations of the given job
		
		subdata=NULL
		for (pop in job_pop){
#			sub=subset(dat,dat[,1]==pop) 08042013 
			if (label==1) sub=subset(dat,dat[,2]==pop)
			if (label==0) sub=subset(dat,dat[,1]==pop)
			subdata=rbind(subdata,sub)
		}
		
		
		
		in_nam=paste('data_job_',id,sep='')
		if (markernames==1) write(mark_line,ncol=length(mark_line),file=in_nam)  # write markernames
		write(t(subdata),ncol=length(subdata[1,]),file=in_nam,append=T)  # write data subset with unique task name-tag
		# generate parameted file for structure 
		param_nam=paste('parameter_job_',as.character(id),sep='')
		out_nam=paste(outpath,'results_job_',as.character(id),sep='')
		
		if (onerowperind==0) nind_job=(length(subdata[,1]))/2
		if (onerowperind==1) nind_job=length(subdata[,1])
		# make list of local parameters 
		
		LocPar=list(name_param=param_nam,outfile=out_nam,infile=in_nam,numinds=nind_job,maxpop=k,burnin=burnin,iter=iter)
		edit_params(GlobPar=GlobPar,LocPar=LocPar)  # edit the parameter file for this particular job
		instr=paste(structure_path,'structure -m ',param_nam,sep='')
		system(instr,ignore.stdout=T)
		
		
		instr=paste('rm ',param_nam,sep='')   # erase parameter file and sub_data file for this job 
		system(instr)
		instr=paste('rm ',in_nam,sep='')
		system(instr)
		
		nr=job$id
		m1=paste('job ',as.character(nr),sep='')
		m2=paste('was calculated ','',sep='')
		tt=Sys.time()
		fff=paste(m1,m2,'at',tt,sep=' ')
		print(fff)
		
				
				
	}
			
			slave_fun_windows=function(job){
#Mac_path="/Applications/Structure.app/Contents/Resources/Java/bin/"	
				structure_path=job$path
				id=job$id
				job_pop=job$job_pop   # list populaiton for this job
				k=as.character(job$k)
				burnin=as.character(job$burnin)
				iter=as.character(job$iter)
				
########################
				temp_nam=paste('temp_list_job_',as.character(id),sep='')
				write(job_pop,file=temp_nam)
				job_pop=as.matrix(read.table(temp_nam,sep=','))
				instr=paste('del ',temp_nam,sep='')
				shell(instr)
				
#### create a data subset for the poulations of the given job
				
				subdata=NULL
				for (pop in job_pop){
					if (label==1) sub=subset(dat,dat[,2]==pop)
					if (label==0) sub=subset(dat,dat[,1]==pop)
					subdata=rbind(subdata,sub)
				}
				
				
				
				in_nam=paste('data_job_',id,sep='')
				if (markernames==1) write(mark_line,ncol=length(mark_line),file=in_nam)  # write markernames
				write(t(subdata),ncol=length(subdata[1,]),file=in_nam,append=T)  # write data subset with unique task name-tag
# generate parameted file for structure 
				param_nam=paste('parameter_job_',as.character(id),sep='')
				out_nam=paste(outpath,'results_job_',as.character(id),sep='')
				
				if (onerowperind==0) nind_job=(length(subdata[,1]))/2
				if (onerowperind==1) nind_job=length(subdata[,1])
# make list of local parameters 
				
				LocPar=list(name_param=param_nam,outfile=out_nam,infile=in_nam,numinds=nind_job,maxpop=k,burnin=burnin,iter=iter)
				
				edit_params(GlobPar=GlobPar,LocPar=LocPar)  # edit the parameter file for this particular job
				instr=paste(structure_path,'structure -m ',param_nam,sep='')
				system(instr,ignore.stdout=T)
				
				
				instr=paste('del ',param_nam,sep='')   # erase parameter file and sub_data file for this job 
				shell(instr)
				instr=paste('del ',in_nam,sep='')
				shell(instr)
				
				nr=job$id
				m1=paste('job ',as.character(nr),sep='')
				m2=paste('was calculated ','',sep='')
				tt=Sys.time()
				fff=paste(m1,m2,'at',tt,sep=' ')
				print(fff)
				
				
				
			}
			
			
			
			
#cl=makeCluster(n_cpu)
			
#		mes=paste('made a ',n_cpu,' cluster', sep='')
#		print(mes)
			
#parLapply(cl,tasks,slave_fun)
			
		if (Sys.info()['sysname']=="Windows"){
			slave_fun=slave_fun_windows
		} else {
			slave_fun=slave_fun_unix
		}
			
			
	
			
		mclapply(tasks,slave_fun,mc.cores=n_cpu)
	
#	stopCluster(cl)
			
			
			m2=paste('was calculated ','',sep='')
			tt=Sys.time()
			fff=paste(m2,'at',tt,sep=' ')
			print(fff)
				
	if (plot_output==1){
		list_id=NULL
		for (i in 1:length(ALL_tasks)){
			job=ALL_tasks[[i]]
			id=job$id
			list_id=c(list_id,id)
		}
#job_pop=job$job_pop   # list populaiton for this job
#			temp_nam=paste('temp_list_job_',as.character(id),sep='')
#			write(job_pop,file=temp_nam)
#			job_pop=as.matrix(read.table(temp_nam,sep=','))
#			instr=paste('rm ',temp_nam,sep='')
#			system(instr)
#			out_nam=paste(outpath,'results_job_',as.character(id),sep='')

		plot_str(loc=outpath,list_id,20,8)
			
	}
}
