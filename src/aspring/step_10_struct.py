import argparse


# get the exon coordinates from the PIR sequence file
def get_sexon_coord(gid,tid,seqFname):

	# open, read and close the input file
	fseq = open(seqFname)
	lines = fseq.readlines()
	fseq.close()
	# go through the file until finding the tid
	i = 0
	found = False
	while (i < len(lines)) and (not found):
		found = lines[i].startswith('>P1;'+gid+' '+tid)
		i = i + 1
	# if the tid was found
	if found:
		sex = lines[i][:-1]
		# should be a list since they are ordered!
		l = []
		sexRef = sex[0]
		startRef = 0
		# go through the whole sequence
		for i in range(1,len(sex)):
			# if the sexon has just changed
			if sex[i] != sexRef:
				# append the previous sexon to the list
				# +1 to the start to shift
				# nothing to the end because we want the one before last
				l.append((sexRef,startRef+1,i))
				sexRef = sex[i]
				startRef = i
		# don't forget the last one!
		# id only one amino acid, startRef+1 should be the last position
		# and i should be equal to len(sex)
		l.append((sexRef,startRef+1,i+1))
		return l


# get the sexon id from the dictionary file
def get_sexon_id(dictFname):

	# open, read and close the input file
	fdic = open(dictFname)
	lines = fdic.readlines()
	fdic.close()
	d = {}
	for line in lines:
		words = line[:-1].split()
		# key: symbol, value: id
		d[words[1]] = words[0]
	print(d)
	return d

def get_asrus(asruFname):
	# res is a dico with as keys the sexons and as values some instance ids
	res = {}
	# readf the data
	fasru = open(asruFname)
	lines = fasru.readlines()
	fasru.close()
	iASRU = 1
	iInst = 1
	# each line should be an asru
	for line in lines[1:]:
		words = line[:-1].strip().split('\"')
		# get the instances
		instances = words[1][1:-1].strip().split(',')
		# for each instance
		for insta in instances:
			# remove extra character
			insta = insta.strip().strip("'")
			# get a list of s-exons
			wds = insta.split(".")
			# if there are more than 1 s-exon
			if len(wds)>1:
				insta = []
				for wd in wds:
					res[wd.strip("$")] = 'a'+str(iASRU)+'i'+str(iInst)
			else:
				res[insta] = 'a'+str(iASRU)+'i'+str(iInst)
			iInst = iInst + 1
		iASRU = iASRU+1
		iInst = 1
	return res

def get_pdb_span(pdb):

	fpdb = open(pdb)
	lines = fpdb.readlines()
	fpdb.close()
	span = []
	for line in lines:
		if line.startswith("ATOM"):
			if line[12:16].strip() == "CA":
				span.append(int(line[22:26]))
	return span

# write out the PML file
def write_pml(coord,d,asrus,pdb,outname):
	fout = open(outname,"w")
	fout.write('load '+pdb+'\n')
	fout.write('bg_color white\n')
	span = get_pdb_span(pdb)
	nameProt = pdb[:-4]
	fout.write('color white '+nameProt+'\n')
	partnerId = 'p'+pdb[-5]
	mycolsex = ['skyblue','yelloworange']
	mycolasru = ['firebrick','lime']
	coli = 0
	colj = 1
	currentInst = ['',[]]
	# go over all sexons
	print(coord)
	for tup in coord:
		start = tup[1]
		end = tup[2]
		# focus on the relevant span
		if start >= span[0] or end <= span[-1]:
			sex = d[tup[0]]
			fout.write('select '+partnerId+'_'+sex+', resid '+str(start)+':'+str(end)+' and '+nameProt+'\n')
			# if sex is in an instance
			if sex in asrus:
				print(sex)
				# different from the current one
				if asrus[sex] != currentInst[0]:
					# if it is not the very first one (current is dummy)
					if currentInst[0] != '':
						# select and color the current one
						namSel = partnerId+currentInst[0]+'_'+'+'.join(currentInst[1])
						fout.write('select '+namSel+', '+' or '.join(currentInst[1])+' and '+nameProt+'\n')
						fout.write('color '+mycolasru[colj]+', '+namSel+'\n')
					# change the current one 
					currentInst[0] = asrus[sex]
					currentInst[1] = [partnerId+'_'+sex]
					# switch instance color
					colj = 1 - colj
				# same as the current one
				else:
					# simply append the sexon to the instance definition
					 currentInst[1].append(partnerId+'_'+sex)
			else:
				# show s-exon color if it's not part of an instance
				fout.write('color '+mycolsex[coli]+', '+partnerId+'_'+d[tup[0]]+'\n')
			# in any case we need to switch sexon color 
			coli = 1 - coli
	# need to select and color the very last one (that won't be replaced anyway...)
	# I'm not sure that is used at all since we will never have something dummy here...
	# maybe if there is only one...?
	if currentInst[0] != '':
		# name of the current instance (s-repeat)
		namSel = partnerId+currentInst[0]+'_'+'+'.join(currentInst[1])
		# select the corresponding s-exon combination
		fout.write('select '+namSel+', '+' or '.join(currentInst[1])+' and '+nameProt+'\n')
		# set the color
		fout.write('color '+mycolasru[colj]+', '+namSel+'\n')
	fout.close()


if __name__ == "__main__":  # execute only if run as a script

    def arg_parser():  # Parser implementation
        parser = argparse.ArgumentParser(prog='highlight_sexons_pml.py',
                                         epilog="   highlight_sexons_pml.py --gidRef ENSG00000029534 --gid ENSG00000029534 --tid ENST00000347528 --pdb ENST00000347528.pdb --asru ENSG00000029534_ASRUs_table.csv",
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('--gid', help='Ensembl gene identifier, eg: ENSG00000029534')        
        parser.add_argument('--tid', help='Ensembl transcript identifier, eg: ENST00000347528')
        parser.add_argument('--gidRef', help='Ensembl reference gene identifier, eg: ENSG00000029534', default='')
        parser.add_argument('--pdb', help='PDB file containing the 3D model of the proteoform')
        parser.add_argument('--asru', help='CSV file containing the ASRUs', default='')
        return parser.parse_args()

    # parse and config arguments
    args = arg_parser()
    if args.gidRef=='':
    	args.gidRef = args.gid
    if args.asru=='':
    	asruFname = args.gidRef+'/'+args.gidRef+'_ASRUs_table.csv'
    else:
    	asruFname = args.asru
    seqFname = args.gidRef+'/thoraxe/phylosofs/transcripts.pir'
    dictFname = args.gidRef+'/thoraxe/phylosofs/s_exons.tsv'
    #print(seqFname,dictFname,args.tid,args.pdb)
    coord = get_sexon_coord(args.gid,args.tid,seqFname)
    d = get_sexon_id(dictFname)
    asrus = get_asrus(asruFname)
    print(asrus)
    write_pml(coord,d,asrus,args.pdb,args.tid+'.pml')