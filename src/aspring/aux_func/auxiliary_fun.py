import os

def get_immediate_subdirectories(a_dir): #while inside directory get subdirectories names into list
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def read_fasta(fasta_file, keep_annotation=False):
    """Read sequences from fasta file.

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.
    keep_annotation : boolean
        Determine is sequence id should contain annotation.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    """

    if not os.path.exists(fasta_file):
        raise InputFileError('Input file %s does not exist.' % fasta_file)

    if os.stat(fasta_file).st_size == 0:
        return {}

    try:

        if fasta_file.endswith('.gz'):
            file_f, file_mode = gzip.open, 'rt'
        else:
            file_f, file_mode = open, 'r'

        seqs = {}
        with file_f(fasta_file, file_mode) as f:

            for line in f.readlines():
                # skip blank lines
                if not line.strip():
                    continue

                if line[0] == '>':
                    if keep_annotation:
                        seq_id = line[1:-1]
                    else:
                        seq_id = line[1:].split(None, 1)[0]

                    seqs[seq_id] = []
                else:
                    seqs[seq_id].append(line.strip())

        for seq_id, seq in seqs.items():
            seqs[seq_id] = ''.join(seq).replace(' ', '')
    except:
        print(traceback.format_exc())
        print()
        print("[Error] Failed to process sequence file: " + fasta_file)
        sys.exit(1)

    return seqs