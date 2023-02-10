# test if seqinr and argparser are installed, if they are not, install them and
# load them

if(!require(seqinr)){
	install.packages("seqinr", dependencies = TRUE, version = "4.2.16")
}

if(!require(argparser)){
	install.packages("argparser", dependencies = TRUE, version = "0.7.1")
}

library("seqinr")
library("argparser")

# Create a parser
p <- arg_parser("get variables given by main")
# Add command line arguments
p <- add_argument(p, "--gene", help="name of gene", type="character")
p <- add_argument(p, "--path_data", help="path to dir containing Thoraxe outputs", type="character")
# Parse the command line arguments
argv <- parse_args(p)
gene_name=argv$gene
path_data = paste0(argv$path_data,"/")

#TO DO : instead read file species_list.txt
species = c("homo_sapiens","gorilla_gorilla","macaca_mulatta","monodelphis_domestica","rattus_norvegicus","mus_musculus","bos_taurus","sus_scrofa","ornithorhynchus_anatinus","xenopus_tropicalis","danio_rerio","caenorhabditis_elegans")
n =  length(species)

combi2<-function(x){
	return(x*(x-1)/2)
}

computeSumOfPairs<-function(fname, match = 1, mismatch = -0.5, gap = 0){
	ali = as.matrix.alignment(read.alignment(fname,format = "FASTA"))
	nspe = dim(ali)[[1]]
	npos = dim(ali)[[2]]
	Ctot = nspe*(nspe-1)/2
	score = 0
	for(j in 1:npos){
		t = table(ali[,j])
		selMatches = which(names(t)!="-"&t>1)
		if(length(selMatches)>0){Cmatch = sum(combi2(t[selMatches]))}
		else{Cmatch = 0}
		if(sum(names(t)=="-")==1){Cgaps = combi2(t["-"]) + t["-"]*(nspe-t["-"])}
		else{Cgaps = 0}
		score = score + Cmatch * match + Cgaps * gap + (Ctot - Cgaps - Cmatch) * mismatch 
	}
	return(score)
}

# get the length of the s-exons in aas, grouping them according to conservation
getSexSizeEvents<-function(){
	genes=c(gene_name)
	#genes = read.table('/Users/antoineszatkownik/Documents/projetAS/human_proteome_1_to_1/benchmark_list.txt',colClass="character")$V1
	path = path_data
	mySizes = c()
	myGenes = c()
	mySpe = c()
	mySpeOcc = c()
	myHum = c()
	mySe = c()
	mySpeTot = c()
	myConsPath = c()
	myTypeEvent = c()
	myScores = c()
	# for each gene
	for(g in genes){
		print(c("treating gene",g))
		# read the s-exon table provided by ThorAxe
		sex = read.table(paste0(path,g,"/thoraxe/s_exon_table.csv"),sep=",",head=TRUE)
		print("s-exon table successfully read")
		# number of species where each s-exon is present
		nSpe = tapply(sex[,"Species"],sex[,"S_exonID"],f<-function(x){return(length(unique(x)))})
		# is each s-exon present in Human?
		isHum = tapply(sex[,"Species"],sex[,"S_exonID"],f<-function(x){return("homo_sapiens"%in%unique(x))})
		# get unique list of s-exon sequences
		sex = unique(sex[,c("Species","S_exonID","S_exon_Start","S_exon_End")])
		# compute their lengths (in aas)
		leSex = sex[,"S_exon_End"]-sex[,"S_exon_Start"]+1
		# compute the max number of residues as the length of the s-exon
		nbaas = tapply(leSex,sex[,"S_exonID"],max)
		# get the number of s-exons
		nse = length(nbaas)
		# order the variables 
		nbaas = nbaas[order(names(nbaas))]
		nSpe = nSpe[order(names(nSpe))]
		isHum = isHum[order(names(isHum))]
		# compute theoretical max MSA scores
		maxScores = nSpe * (nSpe-1) * nbaas / 2
		# add the s-exon sizes to the object mySizes (ordered by the names of the s-exons)
		mySizes = c(mySizes, nbaas)
		# add the gene name to the object myGenes (repeated as many times as there are s-exons in the gene)
		myGenes = c(myGenes, rep(g, nse))
		# add the number of species where each s-exon is present (ordered by the names of the s-exons) 
		mySpe = c(mySpe, nSpe)
		# add the Human booleans (ordered by the names of the s-exons) 
		myHum = c(myHum, isHum) 
		# add the names (ids) of the s-exon in ascending order
		namSe = sort(names(nbaas))
		mySe = c(mySe, namSe)
		# add the number of species where one-to-one orthologs of the gene could be found (repeated as many times as there are s-exons in the gene)
		mySpeTot =c(mySpeTot, rep(length(unique(sex[,"Species"])), nse))
		# get the names of the fasta files
		fnames = paste0(path,g,"/thoraxe/msa/","msa_s_exon_",namSe,".fasta")
		# create a filter for zero s-exons
		#filterZero = startsWith(namSe,"0_")
		# compute sum-of-pair scores of the s-exons that are not zero
		msaScores = rep(0,nse)
		if(sum(maxScores>0)>0){
			msaScores[maxScores>0] = tapply(fnames[maxScores>0], 1:sum(maxScores>0), f<-function(x){mtry = try(read.table(x));if(!inherits(mtry, "try-error")){return(computeSumOfPairs(x))}else{return(0)}}) / maxScores[maxScores>0]}
		myScores = c(myScores, msaScores)
		# read the AS event table provided by ThorAxe
		mtry = try(read.table(paste0(path,g,"/thoraxe/ases_table.csv"),sep=",",head=TRUE,colClass=c(rep("character",6),rep("numeric",2),rep("character",3),rep("numeric",6),rep("character",2))))
		if(!inherits(mtry, "try-error")){
			ases = read.table(paste0(path,g,"/thoraxe/ases_table.csv"),sep=",",head=TRUE,colClass=c(rep("character",6),rep("numeric",2),rep("character",3),rep("numeric",6),rep("character",2)))
			print("AS event table successfully read")
			# get number of events
			nAses = dim(ases)[[1]]}
		else{
			nAses = 0}
		# if there is at least one event
		if(nAses>0){
			# get the number of genes where the alternative path is present
			consAlt = tapply(ases[,"AlternativePathGenes"],1:nAses,f<-function(x){length(strsplit(x,"/")[[1]])})
			# get the number of genes where the canonical path is present
			consCan = tapply(ases[,"CanonicalPathGenes"],1:nAses,f<-function(x){length(strsplit(x,"/")[[1]])})
			# create a data frame 2x3 with the canonical and alternative paths, the respective numbers of species and the type of the ASE (repeated)
			evInfo = data.frame(path=c(ases[,"CanonicalPath"],ases[,"AlternativePath"]),cons=c(consCan,consAlt), aseType = rep(ases[,"ASE"],2))
			# convert event type to character
			evInfo[,"aseType"] = as.character(evInfo[,"aseType"])
			# for each s-exon
			for(i in 1:nse){
				# get its name
				se = namSe[i]
				# add some separator before and after the name
				query  = paste0("/",se,"/")
				# look for the s-exon in the data frame
				if(length(grep(query,evInfo[,1]))>0){
					# get the maximum number of species over all the events involving the s-exon 
					c = max(evInfo[grep(query,evInfo[,1]),"cons"])
					# get the corresponding event type
					d = evInfo[grep(query,evInfo[,1]),"aseType"][which.max(evInfo[grep(query,evInfo[,1]),"cons"])]#toString(sort(unique(ases[grep(query,ases[,1]),"ASE"])))
				}
				# if the s-exon has no event
				else{c = 0;d="no"}
				# add the conservation (number of species) and associated event to the objects myConsPath and myTypeEvent
				myConsPath = c(myConsPath, c)
				myTypeEvent = c(myTypeEvent, d)
			}
		}
		else{
			myConsPath = c(myConsPath, rep(0,nse))
			myTypeEvent = c(myTypeEvent, rep("no",nse))
		}
	}
	# create the final data frame
	res = data.frame(gene = myGenes, speciesTot = mySpeTot,sexon= mySe, size = mySizes, species = mySpe, isHuman = myHum, score = myScores, maxSpeciesASE = myConsPath, aseType = myTypeEvent)
	res[,"size"] = as.numeric(res[,"size"])
	print(c(dim(res)[[1]],sum(res[,"species"]==1),sum(res[,"species"]>1),sum((res[,"size"]<4)&(res[,"species"]>1))))
	#print(res[res[,"size"]<4&res[,"species"]>1,c("gene","sexon")])
	#return(res)
	write.table(res,paste0(path_data,"data/",gene_name,"/",gene_name,"_sexSizeEvents.csv"),sep=",",quote=FALSE,row.names=FALSE)
}

getSexSizeEvents()

# determine which pairs are valid, i.e. not separated by species
checkDuplications<-function(){
	duplications=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_formated.csv"),head=TRUE,sep=",",colClass=c(rep("character",2),rep("numeric",4),"character",rep("numeric",4)))
	duplications=duplications[!duplicated(duplications[,c('S_exon_A','S_exon_B','Gene')]),] #keep only best, so convert one to many table to one to one
	genes = unique(duplications[,"Gene"])
	print(genes)
	path = path_data
	isValid = c()
	for(g in genes){
		#print(c("treating gene",g))
		# read the s-exon table provided by ThorAxe
		sex = read.table(paste0(path,g,"/thoraxe/s_exon_table.csv"),sep=",",head=TRUE)
		#print("s-exon table successfully read")
		subDup= duplications[duplications[,"Gene"]==g,]
		n = dim(subDup)[[1]]
		for(i in 1:n){
			speA = unique(sex[sex[,"S_exonID"]==subDup[i,"S_exon_A"],"Species"])
			speB = unique(sex[sex[,"S_exonID"]==subDup[i,"S_exon_B"],"Species"])
			isValid = c(isValid, length(intersect(speA,speB))>0)
			#if(subDup[i,"S_exon_A"]=='11_2'&subDup[i,"S_exon_B"]=='11_0'){print(g)}
		}
	}
	duplications = cbind(duplications, isValid)
	write.table(duplications,paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid.csv"),sep=",",quote=FALSE,row.names=FALSE)
	return(isValid)
}

checkDuplications()

# determine whether the s-exon is in the canonical path, in the alternative one, or serves as an anchor
determineStatus<-function(sexon, event){
	query  = paste0("/",sexon,"/")
	if(length(grep(query,event["CanonicalPath"]))>0){return("can")}
	if(length(grep(query,event["AlternativePath"]))>0){return("alt")}
	canWords = strsplit(as.character(event["CanonicalPath"]),"/")[[1]]
	anchors = c(canWords[1],canWords[length(canWords)])
	if(sexon%in%anchors){return("anc")}
	return("no")
}

getEventInfo2<-function(){
	duplications=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid.csv"),head=TRUE,sep=",",colClass=c(rep("character",2),rep("numeric",4),"character",rep("numeric",4)))
	genes = unique(duplications[,"Gene"])
	path = path_data
	f<-function(x){return(length(strsplit(x,"/")[[1]])-2)}
	statusDic=c("can","alt")
	res = c()
	for(g in genes){
		print(c("treating gene",g))
		# read the s-exon table provided by ThorAxe
		ases = read.table(paste0(path,g,"/thoraxe/ases_table.csv"),sep=",",head=TRUE,colClass=c(rep("character",6),rep("numeric",2),rep("character",3),rep("numeric",6),rep("character",2)))
        print(ases)
		print("AS event table successfully read")
		subDup = duplications[duplications[,"Gene"]==g,]
		n = dim(subDup)[[1]]
		for(i in 1:n){
			queryA  = paste0("/",subDup[i,"S_exon_A"],"/")
			queryB  = paste0("/",subDup[i,"S_exon_B"],"/")
			indEv = unique(c(grep(queryA,ases[,"CanonicalPath"]),grep(queryA,ases[,"AlternativePath"]),grep(queryB,ases[,"CanonicalPath"]),grep(queryB,ases[,"AlternativePath"])))
            print(paste0('indEv','-',queryA,'-',queryB,'-',indEv))
            print(length(indEv))
			if(length(indEv)>0){
				for(k in indEv){
					tmp = c(g,subDup[i,"S_exon_A"],subDup[i,"S_exon_B"],k,"","","",0,0,FALSE,subDup[i,c("P_value","Cols","Length_A","Length_B")])
					event = ases[k,]
					tmp[5] = event["ASE"]
					# where is the first s-exon ? "no": nowhere, "can": canonical path, "alt": alternative path, "anc": anchor
					tmp[6] = determineStatus(subDup[i,"S_exon_A"],event) #peut être ancre
					tmp[7] = determineStatus(subDup[i,"S_exon_B"],event)
					if(tmp[6]=="can"|tmp[6]=="alt"){
						tmp[8] = length(strsplit(as.character(event[which(statusDic==tmp[6])]),"/")[[1]])-2 #calcul longeur du chemin;-2 car enleve les ancres
					}
					if(tmp[7]=="can"|tmp[7]=="alt"){
						tmp[9] = length(strsplit(as.character(event[which(statusDic==tmp[7])]),"/")[[1]])-2 
					}
					tmp[10] = as.character(event[,"MutualExclusivity"])=="mutually_exclusive"
                    print(tmp)
					res = rbind(res, tmp)
				}
			}
		}
	}
	if(is.null(res)==TRUE){
	  res = data.frame(matrix(nrow = 0, ncol = 14))
	  colnames(res) = c("gene","sexA","sexB","rank","type","statusA","statusB","lePathA","lePathB","exclu","pval","ncols","leA","leB")
	  write.table(res,paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid_analEvents.csv"),sep=",",quote=FALSE,row.names=FALSE)
	}else{
	  colnames(res) = c("gene","sexA","sexB","rank","type","statusA","statusB","lePathA","lePathB","exclu","pval","ncols","leA","leB")
	  write.table(res,paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid_analEvents.csv"),sep=",",quote=FALSE,row.names=FALSE)}

}

getEventInfo2()

treatTrash<-function(){
	eventsDup=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid_analEvents.csv"),head=TRUE,sep=",",colClass=c(rep("character",3),"numeric",rep("character",3),rep("numeric",2),"logical",rep("numeric",4)))
	# the following steps are probably of unnecessary complexity
	selMEHE = which(((eventsDup[,"statusA"]=="alt"&eventsDup[,"statusB"]=="can")|(eventsDup[,"statusA"]=="can"&eventsDup[,"statusB"]=="alt"))&eventsDup[,"exclu"])
	pairsMEHE=apply(unique(eventsDup[selMEHE,1:3]),1,toString)
	filtNotMEHE = !apply(eventsDup[,1:3],1,toString)%in%pairsMEHE
	filtEHE = (eventsDup[,"statusA"]=="alt"&eventsDup[,"statusB"]=="can")|(eventsDup[,"statusA"]=="can"&eventsDup[,"statusB"]=="alt")
	selEHE = which(filtEHE&filtNotMEHE)
	pairsM.EHE=apply(unique(eventsDup[c(selMEHE,selEHE),1:3]),1,toString)
	filtNotM.EHE = !apply(eventsDup[,1:3],1,toString)%in%pairsM.EHE
	# Now I'd like to get those with an anchor
	filtWithAnc = eventsDup[,"statusA"]=="anc"|eventsDup[,"statusB"]=="anc"
	selWithAnc = which(filtNotM.EHE&filtWithAnc)
	pairsM.EHE.Anc = apply(unique(eventsDup[c(selMEHE,selEHE,selWithAnc),1:3]),1,toString)
	filtNotM.EHE.Anc = !apply(eventsDup[,1:3],1,toString)%in%pairsM.EHE.Anc
	filtWithNo = eventsDup[,"statusA"]=="no"|eventsDup[,"statusB"]=="no"
	selWithNo = which(filtNotM.EHE.Anc&filtWithNo)
	for(k in selWithNo){
		words = strsplit(paths[eventsDup[k,"gene"]],"/")[[1]]
		if(eventsDup[k,"statusA"]=="no"){
			if(!eventsDup[k,"sexA"]%in%words){
				filtWithNo[k] = FALSE
			}
		}
		if(eventsDup[k,"statusB"]=="no"){
			if(!eventsDup[k,"sexB"]%in%words){
				filtWithNo[k] = FALSE
			}
		}
	}
	selWithNo = which(filtNotM.EHE.Anc&filtWithNo)
	pairsM.EHE.Anc.No = apply(unique(eventsDup[c(selMEHE,selEHE,selWithAnc,selWithNo),1:3]),1,toString)
	filtNotM.EHE.Anc.No = !apply(eventsDup[,1:3],1,toString)%in%pairsM.EHE.Anc.No
	sel = which(filtNotM.EHE.Anc.No)	
	# unique number of genes displaying these pairs and events
	print(c("unique #(genes) with these pairs and events",length(unique(eventsDup[sel,"gene"]))))
	# unique number of pairs in the set
	print(c("unique #(pairs) in the set",dim(unique(eventsDup[sel,1:3]))[[1]]))
	# unique number of events in the set
	print(c("unique #(events) in the set",dim(unique(eventsDup[sel,c("gene","rank")]))[[1]]))
}



treatFar<-function(withAnc=TRUE){
	sexSizeEvents=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_sexSizeEvents.csv"),sep=",",head=TRUE,row.names=NULL,colClass=c("character","numeric","character",rep("numeric",2),"logical",rep("numeric",2),"character"))
	eventsDup=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid_analEvents.csv"),head=TRUE,sep=",",colClass=c(rep("character",3),"numeric",rep("character",3),rep("numeric",2),"logical",rep("numeric",4)))
	# the following steps are probably of unnecessary complexity
	selMEHE = which(((eventsDup[,"statusA"]=="alt"&eventsDup[,"statusB"]=="can")|(eventsDup[,"statusA"]=="can"&eventsDup[,"statusB"]=="alt"))&eventsDup[,"exclu"])
	pairsMEHE=apply(unique(eventsDup[selMEHE,1:3]),1,toString)
	filtNotMEHE = !apply(eventsDup[,1:3],1,toString)%in%pairsMEHE
	filtEHE = (eventsDup[,"statusA"]=="alt"&eventsDup[,"statusB"]=="can")|(eventsDup[,"statusA"]=="can"&eventsDup[,"statusB"]=="alt")
	selEHE = which(filtEHE&filtNotMEHE)
	pairsM.EHE=apply(unique(eventsDup[c(selMEHE,selEHE),1:3]),1,toString)
	filtNotM.EHE = !apply(eventsDup[,1:3],1,toString)%in%pairsM.EHE
	# Now I'd like to get those with an anchor
	filtWithAnc = eventsDup[,"statusA"]=="anc"|eventsDup[,"statusB"]=="anc"
	selWithAnc = which(filtNotM.EHE&filtWithAnc)
	if(withAnc){
		sel = selWithAnc
		name = "anc"
	}else{
		pairsM.EHE.Anc = apply(unique(eventsDup[c(selMEHE,selEHE,selWithAnc),1:3]),1,toString)
		filtNotM.EHE.Anc = !apply(eventsDup[,1:3],1,toString)%in%pairsM.EHE.Anc
		filtWithNo = eventsDup[,"statusA"]=="no"|eventsDup[,"statusB"]=="no"
		tmp=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_canonical_path.txt"),colClass="character",row.names=1)
		paths=tmp[,1]
		names(paths)=row.names(tmp)
		sel = which(filtNotM.EHE.Anc&filtWithNo)
		print(length(sel))
		for(k in sel){
			words = strsplit(paths[eventsDup[k,"gene"]],"/")[[1]]
			if(eventsDup[k,"statusA"]=="no"){
				if(!eventsDup[k,"sexA"]%in%words){
					filtWithNo[k] = FALSE
				}
			}
			if(eventsDup[k,"statusB"]=="no"){
				if(!eventsDup[k,"sexB"]%in%words){
					filtWithNo[k] = FALSE
				}
			}
		}
		sel = which(filtNotM.EHE.Anc&filtWithNo)
		print(length(sel))
		name = "no"
	}
	# unique number of genes displaying these pairs and events
	print(c("unique #(genes) with these pairs and events",length(unique(eventsDup[sel,"gene"]))))
	# unique number of pairs in the set
	print(c("unique #(pairs) in the set",dim(unique(eventsDup[sel,1:3]))[[1]]))
	# unique number of events in the set
	print(c("unique #(events) in the set",dim(unique(eventsDup[sel,c("gene","rank")]))[[1]]))
	# re-order A and B such that the anchor comes first
	if(identical(sel,integer(0))!=TRUE){
		orderedEvents = t(apply(eventsDup[sel,],1,f<-function(x){indAnc=which(x[c("statusA","statusB")]==name);cbind(x["gene"],x[1+indAnc],x[1+indAnc%%2+1],x[4])}))
		# define three types of s-exon, for each event, anchor, canonical or alternative
		can = rep("",length(sel))
		can[apply(eventsDup[sel,6:7],1,f<-function(x){return("can"%in%x)})] = orderedEvents[apply(eventsDup[sel,6:7],1,f<-function(x){return("can"%in%x)}),3]
		alt = rep("",length(sel))
		alt[apply(eventsDup[sel,6:7],1,f<-function(x){return("alt"%in%x)})] = orderedEvents[apply(eventsDup[sel,6:7],1,f<-function(x){return("alt"%in%x)}),3]
		anc = tapply(orderedEvents[,2],apply(orderedEvents[,c(1,4),drop=FALSE],1,toString),x<-function(x){cf = unique(x);toString(cf[order(nchar(cf), cf)])})
		canNew = tapply(can,apply(orderedEvents[,c(1,4),drop=FALSE],1,toString),x<-function(x){cf = unique(x[x!=""]);toString(cf[order(nchar(cf), cf)])})
		altNew = tapply(alt,apply(orderedEvents[,c(1,4),drop=FALSE],1,toString),x<-function(x){cf = unique(x[x!=""]);toString(cf[order(nchar(cf), cf)])})
		genes = unlist(strsplit(names(anc),", "))[seq(1,2*length(anc),by=2)]
		ranks = unlist(strsplit(names(anc),", "))[seq(2,2*length(anc),by=2)]
		datEvent = data.frame(gene = genes, rank = ranks, anc, can = canNew, alt = altNew)
		print("hello")
		sizeAnc = apply(datEvent,1,f<-function(x){sex=strsplit(x["anc"],", ")[[1]];size=0;for(s in sex){size=size+sexSizeEvents[sexSizeEvents[,"gene"]==x["gene"]&sexSizeEvents[,"sexon"]==s,"size"]};return(size)})
		sizeCan = apply(datEvent,1,f<-function(x){sex=strsplit(x["can"],", ")[[1]];size=0;for(s in sex){size=size+sexSizeEvents[sexSizeEvents[,"gene"]==x["gene"]&sexSizeEvents[,"sexon"]==s,"size"]};return(size)})
		sizeAlt = apply(datEvent,1,f<-function(x){sex=strsplit(x["alt"],", ")[[1]];size=0;for(s in sex){size=size+sexSizeEvents[sexSizeEvents[,"gene"]==x["gene"]&sexSizeEvents[,"sexon"]==s,"size"]};return(size)})
		datEvent = cbind(datEvent, sizeAnc, sizeCan, sizeAlt)
		rownames(datEvent)=1:dim(datEvent)[[1]]
		print("hello")
		ranks = tapply(datEvent[,"rank"],apply(datEvent[,c("gene","anc","can","alt")],1,f<-function(x){paste(x[1],x[2],x[3],x[4],sep=";")}),f<-function(x){if(length(x)>1){return(gsub(" ","",toString(x)))}else{return(x)}})
		datEvent = cbind(unique(datEvent[,-2]),ranks)
		#inAbascal =datEvent[,"gene"]%in%abascal
		#datEvent = cbind(datEvent,inAbascal)
		if(!withAnc){colnames(datEvent)[2] = "no"; colnames(datEvent)[5] = "sizeNo"}
		return(list(sel,datEvent))
	}
	else{}
}


# create a MEHE-specific table indicating the list of mutually exclusive events separating groups of s-exons with some similarity
# the groups of s-exons are given for each event, they may be overlapping, included into one another... across events
# the length of each group is also given
treatEHEs<-function(exclusive=TRUE){
	statusDic=c("can","alt")
	print('a')
	sexSizeEvents=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_sexSizeEvents.csv"),sep=",",head=TRUE,row.names=NULL,colClass=c("character","numeric","character",rep("numeric",2),"logical",rep("numeric",2),"character"))
	print('aa')
	eventsDup=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid_analEvents.csv"),head=TRUE,sep=",",colClass=c(rep("character",3),"numeric",rep("character",3),rep("numeric",2),"logical",rep("numeric",4)))
	# the mutual exckusivity is the most restrictive criterion one can have
	# here we select the pairs and events where one s-exon is in the alternative path, the other in the canonical one
	# and there exists no path in all input transcripts that include both alternative and canonical paths (hence, the 2 s-exons)
	selMEHE = which(((eventsDup[,"statusA"]=="alt"&eventsDup[,"statusB"]=="can")|(eventsDup[,"statusA"]=="can"&eventsDup[,"statusB"]=="alt"))&eventsDup[,"exclu"])
	if(exclusive){
		sel = selMEHE
	}else{
		pairsMEHE=apply(unique(eventsDup[selMEHE,1:3,drop=FALSE]),1,toString)
		filtNotMEHE = !apply(eventsDup[,1:3],1,toString)%in%pairsMEHE
		filtEHE = (eventsDup[,"statusA"]=="alt"&eventsDup[,"statusB"]=="can")|(eventsDup[,"statusA"]=="can"&eventsDup[,"statusB"]=="alt")
		sel = which(filtEHE&filtNotMEHE)
	}
	# unique number of genes displaying these pairs and events
	print(c("unique #(genes) with these pairs and events",length(unique(eventsDup[sel,"gene"]))))
	# unique number of pairs in the set
	print(c("unique #(pairs) in the set",dim(unique(eventsDup[sel,1:3]))[[1]]))
	# unique number of events in the set
	print(c("unique #(events) in the set",dim(unique(eventsDup[sel,c("gene","rank")]))[[1]]))
	# for each event in the selection, one can define the set of s-exon group pairs that are mutually exclusive and share some similarity
	if(identical(sel,integer(0))!=TRUE){
		orderedEvents = t(apply(eventsDup[sel,],1,f<-function(x){indA=which(x["statusA"]==statusDic);cbind(x["gene"],x[1+indA],x[1+indA%%2+1],x[4])}))
		can = tapply(orderedEvents[,2],apply(orderedEvents[,c(1,4),drop=FALSE],1,toString),x<-function(x){cf = unique(x);toString(cf[order(nchar(cf), cf)])})
		print('OOO')
		alt = tapply(orderedEvents[,3],apply(orderedEvents[,c(1,4),drop=FALSE],1,toString),x<-function(x){cf = unique(x);toString(cf[order(nchar(cf), cf)])})
		genes = unlist(strsplit(names(can),", "))[seq(1,2*length(can),by=2)]
		ranks = unlist(strsplit(names(can),", "))[seq(2,2*length(can),by=2)]
		print(ranks)
		datEvent = data.frame(gene = genes, rank = ranks, can = can, alt = alt)
		sizeCan = apply(datEvent,1,f<-function(x){sex=strsplit(x["can"],", ")[[1]];size=0;for(s in sex){size=size+sexSizeEvents[sexSizeEvents[,"gene"]==x["gene"]&sexSizeEvents[,"sexon"]==s,"size"]};return(size)})
		sizeAlt = apply(datEvent,1,f<-function(x){sex=strsplit(x["alt"],", ")[[1]];size=0;for(s in sex){size=size+sexSizeEvents[sexSizeEvents[,"gene"]==x["gene"]&sexSizeEvents[,"sexon"]==s,"size"]};return(size)})
		datEvent = cbind(datEvent, sizeCan, sizeAlt)
		rownames(datEvent)=1:dim(datEvent)[[1]]
		ranks = tapply(datEvent[,"rank"],apply(datEvent[,c("gene","can","alt")],1,f<-function(x){paste(x[1],x[2],x[3],sep=";")}),f<-function(x){if(length(x)>1){return(gsub(" ","",toString(x)))}else{return(x)}})
		datEvent = cbind(unique(datEvent[,-2]),ranks)
		#inAbascal =datEvent[,"gene"]%in%abascal
		#datEvent = cbind(datEvent,inAbascal)
		return(list(sel,datEvent))
	}
	else{}
}
res = treatEHEs()
selMEHE = res[[1]]
datEventMEHE = res[[2]]
res = treatEHEs(FALSE)
selEHE = res[[1]]
datEventEHE = res[[2]]
res = treatFar()
selWithAnc = res[[1]]
datEventWithAnc = res[[2]]
res = treatFar(FALSE)
selWithNo = res[[1]]
datEventWithNo = res[[2]]
eventsDup=read.table(paste0(path_data,"data/",gene_name,"/",gene_name,"_duplication_pairs_filtered_valid_analEvents.csv"),head=TRUE,sep=",",colClass=c(rep("character",3),"numeric",rep("character",3),rep("numeric",2),"logical",rep("numeric",4)))
typePair=rep("No",dim(eventsDup)[[1]])
typePair[selMEHE]="MEX"
typePair[selEHE]="ALT"
typePair[selWithAnc]="REL"
typePair[selWithNo]="UNREL"
write.table(cbind(eventsDup,typePair),paste0(path_data,"data/",gene_name,"/",gene_name,"_eventsDup.txt"),row.names=FALSE,sep=",", quote = FALSE)

# write.table(cbind(rbind(datEventMEHE,datEventEHE),class=c(rep("MEX",dim(datEventMEHE)[[1]]),rep("ALT",dim(datEventEHE)[[1]]))),"../datPairs_MEX_ALT.csv",sep=";",quote=FALSE,row.names=FALSE)
# l1=cbind(datEventWithAnc,rep("REL",dim(datEventWithAnc)[[1]]))
# l2=cbind(datEventWithNo,rep("UNREL",dim(datEventWithNo)[[1]]))
# write.table(l1,"../datPairs_REL.csv",sep=";",quote=FALSE,row.names=FALSE)
# write.table(l2,"../datPairs_UNREL.csv",sep=";",quote=FALSE,row.names=FALSE)
