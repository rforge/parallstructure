
##################### A function to automatize parallel jobs in structure via Rmpi package ############

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


MPI_structure <-
function(joblist=NULL,n_cpu=NULL,structure_path=Mac_path,infile=NULL,outpath=NULL,numinds=NULL,numloci=NULL,get_log=1,
		plot_output=1,label=1,popdata=1,popflag=0,locdata=0,phenotypes=0,markernames=0,mapdist=0,onerowperind=0,phaseinfo=0,
		recessivealleles=0,phased=0,extracol=0,missing=-9,ploidy=2,noadmix=0,linkage=0,usepopinfo=0,locprior=0,
		inferalpha=1,alpha=1.0,popalphas=0,unifprioralpha=1,alphamax=10.0,alphapropsd=0.025,freqscorr=1,onefst=0,
		fpriormean=0.01,fpriorsd=0.05,inferlambda=0,lambda=1.0,computeprob=1,pfromflagonly=0,ancestdist=0,
		startatpopinfo=0,metrofreq=10,updatefreq=1,printqhat=0,revert_convert=0,randomize=1,seed=1234,log10Rstart=NULL,log10Rmin=NULL,log10Rmax=NULL,log10Rpropsd=NULL){
	
	# create a list of global parameters common to all jobs. To be shared with slave nodes
	GlobPar=list(numloci=numloci,label=label,popdata=popdata,popflag=popflag,locdata=locdata,phenotypes=phenotypes,
				 markernames=markernames,mapdist=mapdist,onerowperind=onerowperind,phaseinfo=phaseinfo,recessivealleles=recessivealleles,
				 phased=phased,extracol=extracol,missing=missing,ploidy=ploidy,noadmix=noadmix,linkage=linkage,usepopinfo=usepopinfo,
				 locprior=locprior,inferalpha=inferalpha,alpha=alpha,popalphas=popalphas,unifprioralpha=unifprioralpha,alphamax=alphamax,
				 alphapropsd=alphapropsd,freqscorr=freqscorr,onefst=onefst,fpriormean=fpriormean,fpriorsd=fpriorsd,inferlambda=inferlambda,
				 lambda=lambda,computeprob=computeprob,pfromflagonly=pfromflagonly,ancestdist=ancestdist,startatpopinfo=startatpopinfo,
				 metrofreq=metrofreq,updatefreq=updatefreq,printqhat=printqhat,randomize=randomize,seed=seed,log10Rstart=log10Rstart,log10Rmin=log10Rmin,log10Rmax=log10Rmax,log10Rpropsd=log10Rpropsd)
	if (!is.loaded("mpi_initialize")) {
		library("Rmpi")
	}
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
			if (markernames==1 & recessivealleles==0 & linkage==0){
				mark_line=as.matrix(read.table(data_name,nrows=1))
				recess_line=NULL
				dat=as.matrix(read.table(data_name,skip=1,colClasses="character"))
			}
			
			if (markernames==1 & recessivealleles==1 & linkage==0){
				mark_line=as.matrix(read.table(data_name,nrows=1))
				recess_line=as.matrix(read.table(data_name,skip=1,nrows=1))
				dat=as.matrix(read.table(data_name,skip=2,colClasses="character"))
			}
			
			if (markernames==0 & recessivealleles==0 & linkage==0){
				recess_line=NULL
				mark_line=NULL
				dat=as.matrix(read.table(data_name,colClasses="character"))
			}
			
			if (markernames==0 & recessivealleles==1 & linkage==0){
				recess_line=as.matrix(read.table(data_name,nrows=1))
				mark_line=NULL
				dat=as.matrix(read.table(data_name,skip=1,colClasses="character"))
			}
            
            if (mapdist==1){
                mark_line=as.matrix(read.table(data_name,nrows=1))
                position_line=as.matrix(read.table(data_name,skip=1,nrows=1)) # assumes that if mapdist is given then marker names are given also
                dat=as.matrix(read.table(data_name,skip=2,colClasses="character"))
                print('Warning: in the current implementation, ParallelStructure assumes that if linkage is given, marker names are also given: markernames=1 and recessivealleles=0')
                
            }
            
            if (linkage==1 & mapdist==0) print('Warning: linkage model is specified without mapdistance')
            


#dat=as.matrix(read.table(data_name))
# list_all_pop=unique(dat[,1])  # 08042013
			if (label==1) list_all_pop=unique(dat[,2])
			if (label==0) list_all_pop=unique(dat[,1])
#mark_line=NULL
#if (markernames==1) mark_line=dat[1,]  # import marker names
	Npop=length(list_all_pop)

						
	if (onerowperind==0) nind=(length(dat[,1]))/ploidy
	if (onerowperind==1) nind=length(dat[,1])
							
	
	if (nind!=numinds) {
		m=paste('Error: expected ',numinds,' individuals and found ', nind,' please check data file')
		stop(m)
	}
				
	###############################
						
	#### Import job list ##########
						
	job_list=as.matrix(read.table(joblist,colClasses='character'))
	n_tasks=length(job_list[,1])
	tasks=vector('list')
	for (i in 1:n_tasks){
		job=job_list[i,]
		tasks[[i]]=list(id=as.character(job[1]),job_pop=as.character(job[2]),k=as.character(job[3]),burnin=as.character(job[4]),iter=as.character(job[5]))
		
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
					tasks[[l+1]]=list(id=jobid,job_pop=popi1,k=ki,burnin=burnini,iter=iteri)				}
			}
		}
	}
			
	################################
		ALL_tasks=tasks
	
		
	####### call n_cpu nodes #####
			if (get_log==0) mpi.spawn.Rslaves(nslaves=n_cpu,needlog=FALSE)
			if (get_log==1) mpi.spawn.Rslaves(nslaves=n_cpu)


	#### close nodes in case R crashes ######
	.Last <- function(){
		if (is.loaded("mpi_initialize")){
			if (mpi.comm.size(1) > 0){
				print("Please use mpi.close.Rslaves() to close slaves.")
				mpi.close.Rslaves()
			}
			print("Please use mpi.quit() to quit R")
			.Call("mpi_finalize")
		}
	}


	######### slave function #########


			slave_fun_unix=function(){   #						####### UNIX function
	# Note the use of the tag for sent messages: 
	#     1=ready_for_task, 2=done_task+result, 3=exiting 
	# Note the use of the tag for received messages: 
	#     1=task== go work, 2=no more job-> exit 
		
		
		done='F'
		junk=0
		
		while (done=='F'){
			
			# Signal being ready to receive a new task 
			mpi.send.Robj(junk,0,1) 
			
			# Receive a task 
			job <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
			job_info <- mpi.get.sourcetag() 
			tag <- job_info[2] 
			if (tag==1){  # some job to do
				### gather job infos ###
				#id=as.numeric(job$id)
				id=job$id
				job_pop=job$job_pop   # list populaiton for this job
				k=as.character(job$k)
				burnin=as.character(job$burnin)
				iter=as.character(job$iter)
				
				########################
				#split list of pop_id
				
				T_job_pop=strsplit(job_pop,',')
				job_pop=T_job_pop[[1]]
				
				pop_nr=1:length(job_pop)  
				if (usepopinfo==1){
					# number the populations of the given job from 1 to n as STRUCTURE needs pop id that matches 1..K
					convert=vector('list')
					rev_convert=vector('list')
					for (pop in pop_nr){
						convert[[job_pop[pop]]]=pop	
						rev_convert[[as.character(pop)]]=job_pop[pop]
					}
					
					
				}
				
				
				#### create a data subset for the poulations of the given job
				
				subdata=NULL
				for (pop in job_pop){
#					sub=subset(dat,dat[,1]==pop) 08042013 
					if (label==1) sub=subset(dat,dat[,2]==pop)
					if (label==0) sub=subset(dat,dat[,1]==pop)
					subdata=rbind(subdata,sub)
				}
				
				## replace sub dataset column of popID by converted ids
				if (usepopinfo==1 & label==1){
					ID1=subdata[,2]
					ID2=ID1
					for (i in 1:length(ID1)){
						ID2[i]=convert[[ID1[i]]]
					}
					subdata[,2]=ID2
				}
				if (usepopinfo==1 & label==0){
					ID1=subdata[,1]
					ID2=ID1
					for (i in 1:length(ID1)){
						ID2[i]=convert[[ID1[i]]]
					}
					subdata[,1]=ID2
				}
				############################				
				
				
				## write temprary sub data file
				in_nam=paste('data_job_',id,sep='')
				if (markernames==1 & linkage==0) write(mark_line,ncolumns=length(mark_line),file=in_nam)  # write markernames
				if (recessivealleles==1 & linkage==0) write(recess_line,ncolumns=length(recess_line),file=in_nam,append=T)
                
                if (mapdist==1) {
                    write(mark_line,ncolumns=length(mark_line),file=in_nam)
                    write(position_line,ncolumns=length(position_line),file=in_nam,append=T)
                    
                }
                
				write(t(subdata),ncolumns=length(subdata[1,]),file=in_nam,append=T)  # write data subset with unique task name-tag
				
				# generate parameted file for structure 
				param_nam=paste('parameter_job_',as.character(id),sep='')
				out_nam=paste(outpath,'results_job_',as.character(id),sep='')
				
				if (onerowperind==0) nind_job=(length(subdata[,1]))/GlobPar$ploidy
				if (onerowperind==1) nind_job=length(subdata[,1])
				# make list of local parameters 
				LocPar=list(name_param=param_nam,outfile=out_nam,infile=in_nam,numinds=nind_job,maxpop=k,burnin=burnin,iter=iter)
				
				edit_params(GlobPar=GlobPar,LocPar=LocPar)  # edit the parameter file for this particular job
				instr=paste(structure_path,'structure -m ',param_nam,sep='')
				system(instr)
				
				
				instr=paste('rm ',param_nam,sep='')   # erase parameter file and sub_data file for this job 
				system(instr)
				instr=paste('rm ',in_nam,sep='')
				system(instr)
				
				
				#### reverse the pop ID conversion in output _q and _f file #####
				if (usepopinfo==1 & revert_convert==1){
					ID3=ID1
					if (onerowperind==0){
					
							filter=rep(c(1,0),nind_job)
							ID3=ID1[filter==1]
							
						
					}
					convert_IDs(ID1=ID3,out_nam=out_nam,rev_convert=rev_convert,printqhat=GlobPar$printqhat,k=as.numeric(k))
					
				}
				#################################################################
				out=list(id=id)
				mpi.send.Robj(out,0,2)
				
			}
			else if (tag == 2) { # no more job - goto sleep
				done='T'
				
			}
			

			
		}
		mpi.send.Robj(junk,0,3)	# say goodbye before closing the door

	}

			slave_fun_windows=function(){   #					## WINDOWS function
# Note the use of the tag for sent messages: 
#     1=ready_for_task, 2=done_task+result, 3=exiting 
# Note the use of the tag for received messages: 
#     1=task== go work, 2=no more job-> exit 
				
				
				done='F'
				junk=0
				
				while (done=='F'){
					
					# Signal being ready to receive a new task 
					mpi.send.Robj(junk,0,1) 
					
					# Receive a task 
					job <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
					job_info <- mpi.get.sourcetag() 
					tag <- job_info[2] 
					if (tag==1){  # some job to do
					## gather job infos ###
					#id=as.numeric(job$id)
						id=job$id
						job_pop=job$job_pop   # list populaiton for this job
						k=as.character(job$k)
						burnin=as.character(job$burnin)
						iter=as.character(job$iter)
						
						########################
						#split list of pop_id
						
						T_job_pop=strsplit(job_pop,',')
						job_pop=T_job_pop[[1]]
						
						pop_nr=1:length(job_pop)  
						if (usepopinfo==1){
						# number the populations of the given job from 1 to n as STRUCTURE needs pop id that matches 1..K
							convert=vector('list')
							rev_convert=vector('list')
							for (pop in pop_nr){
								convert[[job_pop[pop]]]=pop
								rev_convert[[as.character(pop)]]=job_pop[pop]
							}
							
							
						}
						
						#### create a data subset for the populations of the given job
						
						subdata=NULL
						
						for (pop in job_pop){
							if (label==1) sub=subset(dat,dat[,2]==pop)
							if (label==0) sub=subset(dat,dat[,1]==pop)
							subdata=rbind(subdata,sub)
						}
						
						## replace sub dataset column of popID by converted ids (only if use popinfo)
						if (usepopinfo==1 & label==1){
							ID1=subdata[,2]
							ID2=ID1
							for (i in 1:length(ID1)){
								ID2[i]=convert[[ID1[i]]]
							}
							subdata[,2]=ID2
						}
						if (usepopinfo==1 & label==0){
							ID1=subdata[,1]
							ID2=ID1
							for (i in 1:length(ID1)){
								ID2[i]=convert[[ID1[i]]]
							}
							subdata[,1]=ID2
						}
						############################	
						
						
						## write temporary sub data file
						in_nam=paste('data_job_',id,sep='')
						if (markernames==1 & linkage==0) write(mark_line,ncolumns=length(mark_line),file=in_nam)  # write markernames
						if (recessivealleles==1 & linkage==0) write(recess_line,ncolumns=length(recess_line),file=in_nam,append=T)
                        
                        if (mapdist==1) {
                            write(mark_line,ncolumns=length(mark_line),file=in_nam)
                            write(position_line,ncolumns=length(position_line),file=in_nam,append=T)
                            
                        }
                        
						write(t(subdata),ncolumns=length(subdata[1,]),file=in_nam,append=T)  # write data subset with unique task name-tag
						# generate parameted file for structure 
						param_nam=paste('parameter_job_',as.character(id),sep='')
						out_nam=paste(outpath,'results_job_',as.character(id),sep='')
						
						if (onerowperind==0) nind_job=(length(subdata[,1]))/GlobPar$ploidy
						if (onerowperind==1) nind_job=length(subdata[,1])
						# make list of local parameters 
						LocPar=list(name_param=param_nam,outfile=out_nam,infile=in_nam,numinds=nind_job,maxpop=k,burnin=burnin,iter=iter)
						
						edit_params(GlobPar=GlobPar,LocPar=LocPar)  # edit the parameter file for this particular job
						instr=paste(structure_path,'structure -m ',param_nam,sep='')
						system(instr)
						
						
						instr=paste('del ',param_nam,sep='')   # erase parameter file and sub_data file for this job 
						shell(instr)
						instr=paste('del ',in_nam,sep='')
						shell(instr)
						
						#### reverse the pop ID conversion in output _q and _f file #####
						if (usepopinfo==1 & revert_convert==1){
							ID3=ID1
							if (onerowperind==0){
								
								filter=rep(c(1,0),nind_job)
								ID3=ID1[filter==1]
								
								
							}
							convert_IDs(ID1=ID3,out_nam=out_nam,rev_convert=rev_convert,printqhat=GlobPar$printqhat,k=as.numeric(k))
							
						}
						#################################################################
						out=list(id=id)
						mpi.send.Robj(out,0,2)
						
					}
					else if (tag == 2) { # no more job - goto sleep
						done='T'
						
					}
					
					
					
				}
				mpi.send.Robj(junk,0,3)	# say goodbye before closing the door
				
			}
			
			
			
			
			
			
			
			
			
			
			
			
	##############################
	############ core program ########
	##############################
	junk=0
	done_slaves=0   # number of unemployed slaves
	n_slaves <- mpi.comm.size()-1
			
	#### share info with slave nodes ######
	mpi.bcast.Robj2slave(GlobPar)
	mpi.bcast.Robj2slave(n_loci)
	mpi.bcast.Robj2slave(outpath)
	mpi.bcast.Robj2slave(onerowperind)
	mpi.bcast.Robj2slave(structure_path)
	mpi.bcast.Robj2slave(dat)
	mpi.bcast.Robj2slave(edit_params)
	mpi.bcast.Robj2slave(markernames)
	mpi.bcast.Robj2slave(mark_line)
	mpi.bcast.Robj2slave(recess_line)
	mpi.bcast.Robj2slave(label)
	mpi.bcast.Robj2slave(usepopinfo)
	mpi.bcast.Robj2slave(printqhat)
	mpi.bcast.Robj2slave(revert_convert)
			
			
	####### Get system info and assign corresponding slave function (Windows/unix)
	if (Sys.info()['sysname']=="Windows"){
		slave_fun=slave_fun_windows
	} else {
		slave_fun=slave_fun_unix
	}
			
	###################################
	mpi.bcast.Robj2slave(slave_fun)
	mpi.bcast.cmd(slave_fun())
	
	##########################


		
			
			
	while (done_slaves < n_slaves){   # until all slaves are done with their jobs
	# Receive a message from a slave 
		message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag()) 
		message_info <- mpi.get.sourcetag() 
		slave_id <- message_info[1] 
		tag <- message_info[2] 
		if (tag == 1) {
			# slave is ready for a task. Give it the next task, or tell it tasks 
			# are done if there are none. 
			if (length(tasks) > 0) {  # there is still job to do
				job=tasks[[1]]
				if (usepopinfo==1){  ####### CHECK concordance of pop_ids with K 
					## pop_id should be 1 2..n; and n ≤ K 
					## otherwise STRUCTURE doesn t take pop info into account :
					##  EG: for K=2 but pop_id=1..3
					## Warning: population prior for individual 201 is 3, which is not
					## in the range 1..2.  Population prior for this individual will
					##	be ignored
					job_pop=job$job_pop   # list populaiton for this job
					k=as.character(job$k)
					T_job_pop=strsplit(job_pop,',')
					job_pop=T_job_pop[[1]]
					if (k<length(job_pop) & popflag==0){
						
						message='Error : if usepopinfo is ==1 the column of prior population id should be 1..n ; and n ≤ K. e.g., 
						if you test K=3 make sure your column of pop ID contains no more than three populations, or use popflag to ignore the prior of some population(s)'
						.Call("mpi_finalize")
						stop(message)
						
					}
					
					
				}
				### if pop id OK, send job to salve node
				Sys.sleep(2)  # wait 2 sec before each job, so the random seed has time to re-new
				mpi.send.Robj(tasks[[1]], slave_id, 1)
	
				tasks[[1]] = NULL
			}
			else {  # no more job to do
				mpi.send.Robj(junk, slave_id, 2) # tell him to close the shop
			}
		}
		else if (tag == 2) { 
	# The message contains results. Do something with the results. 
	# Store them in the data structure
			
			nr=message$id
			m1=paste('job ',as.character(nr),sep='')
			m2=paste('was calculated by node ',as.character(slave_id),sep='')
			tt=Sys.time()
			fff=paste(m1,m2,'at',tt,sep=' ')
			print(fff)
			
		}
		else if (tag == 3) { # a slave have closed down
			done_slaves=done_slaves+1
		}
	}

		
	#print(res)
	mpi.close.Rslaves()
#system('rm *.log')
			
			
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
			
######## generate summary table of results ############
	summary_res=matrix(0,(length(ALL_tasks)+1),8)
	l1=c('run_id','burnin','iterations','k','ln_prob_data','mean_llh','var_llh','mean_alpha')  # column title
	summary_res[1,]=l1
			
			
	for (i in 1:length(ALL_tasks)){
		job=ALL_tasks[[i]]
		res_file=paste(outpath,'results_job_',job$id,'_f',sep='')
		res=as.matrix(read.table(res_file,fill=T,sep='',strip.white=T))
		
		ref_line=NULL
		for (j in 1:length(res[,1])){
			ans=paste(res[j,1],' ',res[j,2],' ',res[j,3],' ',res[j,4],sep='')  # fish for the relevant line in the result file
			if (ans=='Estimated Ln Prob of') ref_line=j
		}
		Ln_P_data=res[ref_line,7]
		mean_llh=res[(ref_line+1),7]
		Var_llh=res[(ref_line+2),6]
		mean_alpha=res[(ref_line+3),6]
		
		
		l=c(job$id,job$burnin,job$iter,job$k,Ln_P_data,mean_llh,Var_llh,mean_alpha)
		summary_res[(i+1),]=l
		
		
	}
	
	write.csv(summary_res,file='results_summary.csv')
######################################################
}
