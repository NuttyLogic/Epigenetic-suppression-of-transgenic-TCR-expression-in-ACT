from BSBolt.Align.AlignmentHelpers import convert_alpha_numeric_cigar


class ProcessVectorSpanningReads:
    """Indentify high quality reads or read pairs that span a vector of interest and the genome."""

    def __init__(self, read_uniqueness_threshold: float = 1.1, read_mapping_threshold: float = 0.9,
                 multibase_threshold: float = 0.1):
        self.first_read = {'65', '67', '73', '81', '89', '97', '113', '115', '321', '323', '345', '329', '369', '371'}
        self.fr_reference = {'W_C2T', 'C_G2A'}
        self.read_uniqueness_threshold = read_uniqueness_threshold
        self.read_mapping_threshold = read_mapping_threshold
        self.multibase_threshold = multibase_threshold

    def get_integration_sites(self, read_group: list = None, vector: str = None):
        first_reads, second_reads = self.get_paired_reads(read_group)
        group_1 = self.process_read_mapping(first_reads)
        group_2 = self.process_read_mapping(second_reads)
        g1_vector, g1_genome, g1_split = self.assess_alignment_contigs(group_1, vector)
        g2_vector, g2_genome, g2_split = self.assess_alignment_contigs(group_2, vector)
        # discard read groups with conflicting integration site information
        if g1_split and g2_split:
            return None
        # discard read group with only chromosome or vector mapping information
        if not any((g1_vector, g2_vector)) or not any((g1_genome, g2_genome)):
            return None
        # discordant reads
        if not g1_split and not g2_split:
            # read information -> QNAME, CHROM, mapping_ref, l_reference_pos, r_reference_pos,
            # l_query_pos, r_query_pos, Alignment_Score, list[matched_base_position: int]
            return self.process_discordant_int(group_1, group_2, g1_genome)
        return self.process_split_int(group_1, group_2, g1_split, g1_genome, g2_genome, vector)

    @staticmethod
    def process_split_int(group_1: list, group_2: list, g1_split: bool, g1_genome: bool, g2_genome: bool, vector: str):
        split_group = group_1 if g1_split else group_2
        supporting_group = group_1 if not g1_split else group_2
        supporting_genome = g1_genome if not g1_split else g2_genome
        genome_split, vector_split = None, None
        for read in split_group:
            if read[1] == vector:
                vector_split = read
            else:
                genome_split = read
        ref_pos = genome_split[3] if genome_split[5] > vector_split[5] else genome_split[4]
        if not supporting_group:
            return 'split_single', genome_split[0], genome_split[1], ref_pos, genome_split[7], vector_split[7]
        else:
            if supporting_genome:
                if supporting_group[0][1] != genome_split[0][1]:
                    return None
                return 'split_paired', genome_split[0], genome_split[1], ref_pos, \
                       genome_split[7] + supporting_group[0][7], vector_split[7]
            else:
                return 'split_paired', genome_split[0], genome_split[1], ref_pos, \
                       genome_split[7], vector_split[7] + supporting_group[0][7]

    def process_discordant_int(self, group_1, group_2, g1_genome):
        read_1 = group_1[0]
        read_2 = group_2[0]
        if g1_genome:
            ref_pos = read_1[4] if read_1[2] in self.fr_reference else read_1[3]
            return 'discord_1', read_1[0], read_1[1], ref_pos, read_1[7], read_2[7]
        else:
            ref_pos = read_2[3] if read_2[2] in self.fr_reference else read_2[4]
            return 'discord_2', read_2[0], read_2[1], ref_pos, read_2[7], read_1[7]

    @staticmethod
    def assess_alignment_contigs(read_group: list, vector: str) -> (bool, bool):
        vector_mapping, genome_mapping = False, False
        genome_contigs = []
        for read in read_group:
            if read[1] == vector:
                vector_mapping = True
            else:
                genome_contigs.append(read[1])
                genome_mapping = True
        if len(genome_contigs) > 1:
            genome_mapping = False
            vector_mapping = False
        return vector_mapping, genome_mapping, all((vector_mapping, genome_mapping))

    def get_paired_reads(self, read_group):
        first_reads, second_reads = [], []
        for read in read_group:
            if read[1] in self.first_read:
                first_reads.append(read)
            else:
                second_reads.append(read)
        return first_reads, second_reads

    def process_read_mapping(self, read_paired):
        processed_reads = []
        for read in read_paired:
            cigar_tuple = convert_alpha_numeric_cigar(read[4])
            left_most_pos = int(read[3])
            # read information -> l_query_pos, r_query_pos, l_reference_pos,
            # r_reference_pos, list[matched_base_position: int]
            l_pos, r_pos, reference_pos, matched_base_pos = self.get_mapped_bases(cigar_tuple, left_most_pos)
            # read information -> QNAME, CHROM, mapping_ref, l_ref_pos, r_ref_pos, l_query_pos,
            # r_query_pos, Alignment_Score, list[matched_base_position: int]
            processed_reads.append(
                (read[0], read[2], read[8], left_most_pos, reference_pos, l_pos, r_pos, read[7], matched_base_pos))
        processed_reads.sort(key=lambda x: x[7], reverse=True)
        if not processed_reads:
            return processed_reads
        if len(processed_reads) < 2:
            return [read[0:8] for read in processed_reads]
        else:
            primary_reads = [processed_reads[0]]
            duplicate_read = False
            for read in processed_reads[1:]:
                for p_read in primary_reads:
                    if self.get_duplication_proportion(read[8], p_read[8]) > self.multibase_threshold:
                        if read[7] == p_read[7]:
                            return []
                        duplicate_read = True
                        break
                if not duplicate_read:
                    primary_reads.append(read)
        return [read[0:8] for read in primary_reads]

    @staticmethod
    def get_duplication_proportion(read_bases, comparison_bases):
        return len([base for base in read_bases if base in comparison_bases]) / len(read_bases)

    @staticmethod
    def get_mapped_bases(cigar_tuple: tuple, reference_position: int) -> (int, int, list):
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
            elif cigar_type not in query_consumers and cigar_type in reference_consumers:
                reference_position += 1
        return left_mapped_pos, right_mapped_pos, reference_position, matched_base_positions
