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
        alignment_score = None
        mapping_reference = None
        for tag in SAM_TAGS:
            if tag[0:3] == 'AS:':
                alignment_score = tag.split(':')[-1]
            elif tag[0:4] == 'XO:Z':
                mapping_reference = tag.split(':')[-1]
        yield QNAME, FLAG, RNAME, RNEXT, int(POS), CIGAR, int(alignment_score), mapping_reference


def get_spanning_reads(file_path: str = None, plasmid_names: set = None) -> dict:
    mapped_reads = {}
    for sam_read in stream_mapped_reads(file_path, included_flag=4, excluded_flag=1540):
        QNAME, FLAG, RNAME, RNEXT, POS, CIGAR, alignment_score, mapping_reference = sam_read
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
            for read in read_group[0]:
                if read[2][0:3] == 'chr':
                    plasmid_reads[qname] = read_group[0]
                    break
    return plasmid_reads


def call_read_integrations(integration_processor, vector_mapped_reads, vector):
    integration_calls = []
    for read_label, read_group in vector_mapped_reads.items():
        called_integration = integration_processor.get_integration_sites(read_group, vector=vector)
        if called_integration:
            integration_calls.append(called_integration)
    return integration_calls
