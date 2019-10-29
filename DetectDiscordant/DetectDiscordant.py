from BSBolt.Align.AlignmentHelpers import convert_alpha_numeric_cigar


class ProcessVectorSpanningReads:
    """Indentify high quality reads or read pairs that span a vector of interest and the genome."""

    def __init__(self, read_uniqueness_threshold: float = 1.1, read_mapping_threshold: float = 0.9,
                 multibase_threshold: float = 0.1):
        self.first_read = {'65', '67', '73', '81', '89', '97', '113', '115', '321', '323', '345', '329', '369', '371'}
        self.proper_pair = {'67', '131', '115', '179', '323', '371', '387', '435'}
        self.read_uniqueness_threshold = read_uniqueness_threshold
        self.read_mapping_threshold = read_mapping_threshold
        self.multibase_threshold = multibase_threshold

    def get_integration_sites(self, read_group: list = None, vector: str = None):
        first_reads, second_reads = self.get_paired_reads(read_group)
        group_1 = self.process_read_mapping(first_reads, vector)
        group_2 = self.process_read_mapping(second_reads, vector)
        read_1_pass, read1_integration_site = self.assess_alignments(group_1, vector)
        read_2_pass, read2_integration_site = self.assess_alignments(group_2, vector)
        print(read_1_pass, read_2_pass)
        print(read1_integration_site, read2_integration_site)

    def process_read_groups(self, group_1, group_2):
        if group_1[6] and group_2[6]:
            return None
        pass

    def assess_alignments(self, group, vector):
        total_bases, multi_bases, vector_bases, genome_bases, read_len, mapped_contigs, vector_mapping = group
        read_uniqueness = (total_bases / read_len) < self.read_uniqueness_threshold
        read_mapped_prop = ((total_bases - multi_bases) / read_len) > self.read_mapping_threshold
        multi_mapping_bases = (multi_bases / read_len) < self.multibase_threshold
        mapping_count = 0
        for mappings in mapped_contigs.values():
            mapping_count += len(mappings)
        contig_mappings = mapping_count == 1
        if vector_mapping:
            contig_mappings = mapping_count <= 2
        if not all((read_uniqueness, read_mapped_prop, multi_mapping_bases, contig_mappings)):
            return None
        integration_site = None
        if vector_mapping and len(mapped_contigs) > 1:
            integration_site = self.get_split_site(mapped_contigs, vector)
        return True, integration_site

    def get_split_site(self, mapped_contigs, vector):
        chrom_site, vector_site = None, None
        integration_contig = None
        chrom_split, vector_split = None, None
        for contig, sites in mapped_contigs.items():
            if contig == vector:
                vector_site = sites[0]
            else:
                integration_contig = contig
                chrom_site = sites[0]
        if not chrom_site:
            return None
        if chrom_site[2] < vector_site[2]:
            chrom_split = chrom_site[1]
            vector_split = vector_site[0]
        else:
            chrom_split = chrom_site[0]
            vector_split = vector_site[1]
        return integration_contig, chrom_split, vector_split

    def get_paired_reads(self, read_group):
        first_reads, second_reads = [], []
        for read in read_group:
            if read[1] in self.first_read:
                first_reads.append(read)
            else:
                second_reads.append(read)
        return first_reads, second_reads

    def process_read_mapping(self, read_paired, vector):
        observed_bases = []
        mapped_contigs = {}
        vector_bases, genome_bases = 0, 0
        read_len = None
        vector_mapping = False
        for read in read_paired:
            cigar_tuple = convert_alpha_numeric_cigar(read[4])
            left_most_pos = int(read[3])
            read_bases, read_len, right_most_pos, l_pos, r_pos = self.get_mapped_bases(cigar_tuple, left_most_pos)
            if read[2] not in mapped_contigs:
                mapped_contigs[read[2]] = [(left_most_pos, right_most_pos, l_pos, r_pos)]
            else:
                mapped_contigs[read[2]].append((left_most_pos, right_most_pos, l_pos, r_pos))
            observed_bases.extend(read_bases)
            if read[2] == vector:
                vector_bases += len(read_bases)
                vector_mapping = True
            else:
                genome_bases += len(read_bases)
        total_bases = len(observed_bases)
        multi_bases = total_bases - len(set(observed_bases))
        return total_bases, multi_bases, vector_bases, genome_bases, read_len, mapped_contigs, vector_mapping

    def get_mapped_bases(self, cigar_tuple, reference_position):
        matched_base_positions = []
        reference_consumers = {0, 2, 3, 7, 8}
        query_consumers = {0, 1, 4, 7, 8}
        # set relative to genomic position so add reference start and one since capturing the first base
        query_position = 0
        left_mapped_pos, right_mapped_pos = None, None
        for cigar_type, cigar_count in cigar_tuple:
            if cigar_type in reference_consumers and cigar_type in query_consumers:
                for _ in range(cigar_count):
                    if not left_mapped_pos:
                        left_mapped_pos = query_position
                    right_mapped_pos = query_position
                    matched_base_positions.append(query_position)
                    query_position += 1
                    reference_position += 1
            elif cigar_type in query_consumers and cigar_type not in reference_consumers:
                query_position += cigar_count
            elif not cigar_type in query_consumers and cigar_type in reference_consumers:
                reference_position += 1
        return matched_base_positions, query_position, reference_position, left_mapped_pos, right_mapped_pos
