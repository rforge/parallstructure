edit_params <-
function(LocPar,GlobPar){


	out=paste(getwd(),'/',LocPar$outfile,sep='')
	l=paste('#define OUTFILE ',out,sep='')
	write(l,file=LocPar$name_param)
	
	infi=paste(getwd(),'/',LocPar$infile, sep='')
	l=paste('#define INFILE ',infi,sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define NUMINDS ',as.character(LocPar$numinds),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define NUMLOCI ',as.character(GlobPar$numloci),sep='')
	write(l,file=LocPar$name_param,append=T)

	l=paste('#define LABEL ',as.character(GlobPar$label),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define POPDATA ',as.character(GlobPar$popdata),sep='')
	write(l,file=LocPar$name_param,append=T)

	l=paste('#define POPFLAG ',as.character(GlobPar$popflag),sep='')
	write(l,file=LocPar$name_param,append=T)

	l=paste('#define LOCDATA ',as.character(GlobPar$locdata),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define PHENOTYPE ',as.character(GlobPar$phenotype),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define MARKERNAMES ',as.character(GlobPar$markernames),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define MAPDISTANCES ',as.character(GlobPar$mapdistances),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define ONEROWPERIND ',as.character(GlobPar$onerowperind),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define PHASEINFO ',as.character(GlobPar$phaseinfo),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define PHASED ',as.character(GlobPar$phased),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define RECESSIVEALLELES ',as.character(GlobPar$recessivealleles),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define EXTRACOLS ',as.character(GlobPar$extracol),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define MISSING ',as.character(GlobPar$missing),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define PLOIDY ',as.character(GlobPar$ploidy),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define MAXPOPS ',as.character(LocPar$maxpop),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define BURNIN ',as.character(LocPar$burnin),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define NUMREPS ',as.character(LocPar$iter),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define NOADMIX ',as.character(GlobPar$noadmix),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define LINKAGE ',as.character(GlobPar$linkage),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define USEPOPINFO ',as.character(GlobPar$usepopinfo),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define LOCPRIOR ',as.character(GlobPar$locprior),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define INFERALPHA ',as.character(GlobPar$inferalpha),sep='')
	write(l,file=LocPar$name_param,append=T)

	
	l=paste('#define ALPHA ',as.character(GlobPar$alpha),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define POPALPHAS ',as.character(GlobPar$popalphas),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define UNIFPRIORALPHA ',as.character(GlobPar$unifprioralpha),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define ALPHAMAX ',as.character(GlobPar$alphamax),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define ALPHAPROPSD ',as.character(GlobPar$alphapropsd),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define FREQSCORR ',as.character(GlobPar$freqscorr),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define ONEFST ',as.character(GlobPar$onefst),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define FPRIORMEAN ',as.character(GlobPar$fpriormean),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define FPRIORSD ',as.character(GlobPar$fpriorsd),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define INFERLAMBDA ',as.character(GlobPar$inferlambda),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define LAMBDA ',as.character(GlobPar$lambda),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define COMPUTEPROB ',as.character(GlobPar$computeprob),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define PFROMPOPFLAGONLY ',as.character(GlobPar$pfromflagonly),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define ANCESTDIST ',as.character(GlobPar$ancestdist),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define STARTATPOPINFO ',as.character(GlobPar$startatpopinfo),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define METROFREQ ',as.character(GlobPar$metrofreq),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define UPDATEFREQ ',as.character(GlobPar$updatefreq),sep='')
	write(l,file=LocPar$name_param,append=T)
	
	l=paste('#define PRINTQHAT ',as.character(GlobPar$printqhat),sep='')
	write(l,file=LocPar$name_param,append=T)

	l=paste('#define RANDOMIZE ',as.character(GlobPar$randomize),sep='')
	write(l,file=LocPar$name_param,append=T)
		
	

	
	
	

}
