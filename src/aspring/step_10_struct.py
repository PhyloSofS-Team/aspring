import shutil
import pandas as pd
import requests
from biomart import BiomartServer
import os
import argparse
import sys


def parse_args(args):
    parser = get_arg_parser()
    return parser.parse_args(args)


def get_arg_parser():
    parser = argparse.ArgumentParser(
        description='STEP 10 : Generates PyMOL scripts to visualize protein structures with highlighted s-exons and ASRUs.'
    )
    parser.add_argument('--gene',
                        dest='geneName',
                        type=str,
                        required=True,
                        help='name of gene')
    parser.add_argument('--dataPATH',
                        type=str,
                        required=True,
                        help='path to dir containing Thoraxe outputs')
    return parser

# get the exon coordinates from the PIR sequence file


def get_sexon_coord(gid, tid, seqFname):

    # open, read and close the input file
    fseq = open(seqFname)
    lines = fseq.readlines()
    fseq.close()
    # go through the file until finding the tid
    i = 0
    found = False
    while (i < len(lines)) and (not found):
        found = lines[i].startswith('>P1;'+gid+' ') and tid in lines[i]
        i = i + 1
    # if the tid was found
    if found:
        sex = lines[i][:-1]
        # should be a list since they are ordered!
        l = []
        sexRef = sex[0]
        startRef = 0
        # go through the whole sequence
        for i in range(1, len(sex)):
            # if the sexon has just changed
            if sex[i] != sexRef:
                # append the previous sexon to the list
                # +1 to the start to shift
                # nothing to the end because we want the one before last
                l.append((sexRef, startRef+1, i))
                sexRef = sex[i]
                startRef = i
        # don't forget the last one!
        # id only one amino acid, startRef+1 should be the last position
        # and i should be equal to len(sex)
        l.append((sexRef, startRef+1, i+1))
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
            if len(wds) > 1:
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

# Function to get the UniProt ID from the AlphaFold DB PDB file name
# For example: A0A2C9C333 from AF-A0A2C9C333-F1-model_v4.pdb


def get_uniprot_id_from_filename(pdb_filename):
    # Split the file name at the hyphen and underscore characters
    parts = pdb_filename.split("-")
    uniprot_id = parts[1]  # Extract the UniProt ID
    return uniprot_id


def write_pml(coord, d, asrus, pdb, outname):
    fout = open(outname, "w")
    fout.write('load '+pdb+'\n')
    fout.write('bg_color white\n')
    span = get_pdb_span(pdb)
    nameProt = pdb.split('/')[-1][:-4]
    fout.write('color white '+nameProt+'\n')
    # partnerId = 'p'+pdb[-5]
    partnerId = get_uniprot_id_from_filename(nameProt)
    mycolsex = ['skyblue', 'yelloworange']
    mycolasru = ['firebrick', 'lime']
    coli = 0
    colj = 1
    currentInst = ['', []]
    # go over all sexons
    print(coord)
    for tup in coord:
        start = tup[1]
        end = tup[2]
        # focus on the relevant span
        if start >= span[0] or end <= span[-1]:
            sex = d[tup[0]]
            fout.write('select '+partnerId+'_'+sex+', resid ' +
                       str(start)+':'+str(end)+' and '+nameProt+'\n')
            # if sex is in an instance
            if sex in asrus:
                print(sex)
                # different from the current one
                if asrus[sex] != currentInst[0]:
                    # if it is not the very first one (current is dummy)
                    if currentInst[0] != '':
                        # select and color the current one
                        namSel = partnerId+currentInst[0]+'_'+'+'.join(currentInst[1])
                        fout.write('select '+namSel+', ' +
                                   ' or '.join(currentInst[1])+' and '+nameProt+'\n')
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
        fout.write('select '+namSel+', ' +
                   ' or '.join(currentInst[1])+' and '+nameProt+'\n')
        # set the color
        fout.write('color '+mycolasru[colj]+', '+namSel+'\n')
    fout.write('deselect\n')
    fout.close()


def download_alphafold_structure(uniprot_id, output_path='.'):
    base_url = 'https://alphafold.ebi.ac.uk/files/'
    file_name = f'AF-{uniprot_id}-F1-model_v4.pdb'
    url = f'{base_url}{file_name}'

    response = requests.get(url)

    output_file_path = os.path.join(output_path, file_name)

    if response.status_code == 200:
        with open(output_file_path, 'wb') as f:
            f.write(response.content)
        print(f'Successfully downloaded {file_name} to {output_file_path}')
    else:
        output_file_path = None
        print(
            f'Error: Could not download {file_name}. HTTP status code: {response.status_code}')

    return output_file_path


def get_transcripts(s_exon_table):
    # Read the CSV file
    df = pd.read_csv(s_exon_table)

    # Initialize an empty list to store the tuples
    result = []

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        # Get the species and modify it as required
        species_full = row['Species'].split('_')
        species = species_full[0][0] + species_full[1]

        # Get the gene ID
        gene_id = row['GeneID']

        # Get the transcript ID cluster and split it by '/'
        transcript_id_cluster = row['TranscriptIDCluster'].split('/')

        # Iterate through the transcript IDs
        for transcript_id in transcript_id_cluster:
            # Create a tuple and append it to the result list
            result.append((species, gene_id, transcript_id))

    # return only the unique tuples
    return list(set(result))


def get_uniprot_ids(species_gene_transcript_ids):
    # Create a DataFrame from the output of process_csv
    df = pd.DataFrame(species_gene_transcript_ids, columns=[
                      'species', 'gene_id', 'transcript_id'])

    # Split transcripts with '/' into different rows
    df = df.assign(transcript_id=df['transcript_id'].str.split(
        '/')).explode('transcript_id')

    # Keep only unique rows
    df = df.drop_duplicates()

    # Group by species
    grouped = df.groupby('species')

    server = BiomartServer("http://www.ensembl.org/biomart")
    uniprot_ids = []

    data = []
    for species, group in grouped:
        dataset = server.datasets[f'{species}_gene_ensembl']

        gene_ids = group['gene_id'].tolist()
        transcript_ids = group['transcript_id'].tolist()

        response = dataset.search({
            'attributes': ['ensembl_gene_id', 'ensembl_transcript_id', 'uniprotswissprot', 'uniprotsptrembl'],
            'filters': {'ensembl_gene_id': gene_ids, 'ensembl_transcript_id': transcript_ids}
        })

        for line in response.iter_lines():
            line = line.decode('utf-8')
            fields = line.split('\t')
            fetched_gene_id, fetched_transcript_id, uniprotswissprot, uniprotsptrembl = fields
            uniprot_ids = [up for up in [uniprotswissprot, uniprotsptrembl]]
            for uniprot_id in uniprot_ids:
                if fetched_gene_id in gene_ids and fetched_transcript_id in transcript_ids:
                    data.append((fetched_gene_id, fetched_transcript_id, uniprot_id))

    return pd.DataFrame(data, columns=['gene_id', 'transcript_id', 'uniprot_id'])


def download_pdb_structures(s_exon_table, output_path='.'):
    structures = []
    species_gene_transcript_ids = get_transcripts(s_exon_table)
    df = get_uniprot_ids(species_gene_transcript_ids)
    for row in df.itertuples():
        gene_id = row.gene_id
        transcript_id = row.transcript_id
        uniprot_id = row.uniprot_id
        try:
            if uniprot_id:
                pdb_file = download_alphafold_structure(
                    uniprot_id, output_path=output_path)
                if pdb_file is not None:
                    structures.append((gene_id, transcript_id, pdb_file))
            else:
                print(
                    f'Error: Could not find UniProt ID for gene ID {gene_id} and transcript ID {transcript_id}')
        except Exception as e:
            print(f'Error: {e}')
    return structures


def get_pdb_sequence_length(pdb):
    sequence_length = 0
    with open(pdb, 'r') as f:
        for line in f.readlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                sequence_length += 1
    return sequence_length


def get_transcript_sequence_length(gid, tid, seqFname):
    with open(seqFname, 'r') as fseq:
        lines = fseq.readlines()

    i = 0
    found = False
    while i < len(lines) and not found:
        found = lines[i].startswith('>P1;' + gid + ' ') and tid in lines[i]
        i += 1

    if found:
        return len(lines[i].strip())
    else:
        return None


def run():
    args = parse_args(sys.argv[1:])
    gene = args.geneName
    path_data = args.dataPATH

    seqFname = os.path.join(path_data, gene, 'thoraxe', 'phylosofs', 'transcripts.pir')
    dictFname = os.path.join(path_data, gene, 'thoraxe', 'phylosofs', 's_exons.tsv')
    s_exon_table = os.path.join(path_data, gene, 'thoraxe', 's_exon_table.csv')
    asruFname = os.path.join(path_data, 'data', gene, f'{gene}_ASRUs_table.csv')

    if os.path.exists(asruFname):
        output_path = os.path.join(path_data, 'data', gene, 'structures')
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.makedirs(output_path)

        structures = download_pdb_structures(s_exon_table, output_path=output_path)

        for gene_id, transcript_id, pdb_file in structures:
            coord = get_sexon_coord(gene_id, transcript_id, seqFname)
            if coord is None:
                print(f'Error: Could not find {transcript_id} (gene ID: {gene_id})')
                continue

            pdb_seq_length = get_pdb_sequence_length(pdb_file)
            transcript_seq_length = get_transcript_sequence_length(
                gene_id, transcript_id, seqFname)
            if pdb_seq_length != transcript_seq_length:
                continue

            d = get_sexon_id(dictFname)
            asrus = get_asrus(asruFname)
            pdb_code = pdb_file.split("/")[-1].split(".")[0]
            script_file = os.path.join(
                output_path, f'{gene_id}_{transcript_id}_{pdb_code}.pml')
            write_pml(coord, d, asrus, pdb_file, script_file)
    else:
        print(f'Error: Could not find {asruFname}')


if __name__ == '__main__':
    run()
