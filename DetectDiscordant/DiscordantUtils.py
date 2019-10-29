import subprocess


def stream_mapped_reads(file_path, included_flag=None, excluded_flag=None):
    """ process reads streamed using samtools view, samtools must be on
        path for this to work """
    stream_command = ['samtools', 'view']
    if included_flag:
        stream_command.extend(['-F', str(included_flag)])
    if excluded_flag:
        stream_command.extend(['-f', str(excluded_flag)])
    stream_command.append(file_path)
    read_stream = subprocess.Popen(stream_command, stdout=subprocess.PIPE, universal_newlines=True)
    for line in iter(read_stream.stdout.readline, ''):
        QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, *SAM_TAGS = line.strip().split('\t')
        yield QNAME, FLAG, RNAME, POS, CIGAR, RNEXT, PNEXT


def get_plasmid_reads(file_path: str = None, plasmid_names: set = None):
    mapped_reads = {}
    for sam_read in stream_mapped_reads(file_path, included_flag=4):
        QNAME, FLAG, RNAME, POS, CIGAR, RNEXT, PNEXT = sam_read
        plasmid_read = RNAME in plasmid_names
        if QNAME not in mapped_reads:
            mapped_reads[QNAME] = [[sam_read], plasmid_read]
        else:
            if plasmid_read:
                mapped_reads[QNAME][0].append(sam_read)
                mapped_reads[QNAME][1] = plasmid_read
            else:
                mapped_reads[QNAME][0].append(sam_read)
    plasmid_reads = {}
    for qname, read_group in mapped_reads.items():
        if read_group[1]:
            plasmid_reads[qname] = read_group
    return plasmid_reads


def return_plasmid_reads(file_path: str = None, plasmid_names: set = None, return_dict: dict = None, sample_name: str = None):
    plasmid_reads = get_plasmid_reads(file_path = file_path, plasmid_names = plasmid_names)
    return_dict[sample_name] = plasmid_reads


def propogate_error(error):
    raise error


def check_cigar_tuple_end(cigar_tuple):
    cigar_pass = {0, 4}
    cigar_start = cigar_tuple[0][0]
    cigar_end = cigar_tuple[-1][0]
    if cigar_start in cigar_pass and cigar_end in cigar_pass:
        if cigar_start != cigar_end:
            return True
    return False

