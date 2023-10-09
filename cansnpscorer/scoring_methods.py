from abc import ABC, abstractmethod

class ScoringMethod(ABC):
    def __init__(self, name):
        self.name = name

    @abstractmethod
    def calculate(self, snp, snps, reference, genes, tandem_repeats):
        pass

class ScoringManager:
    def __init__(self, methods):
        self.methods = methods

    def calculate_scores(self, snps, reference, genes, tandem_repeats):
        scores = {}
        for snp in snps:
            snp_score = 0
            for method in self.methods:
                snp_score += method.calculate(snp, snps, reference, genes, tandem_repeats)
            scores[snp] = snp_score
        return scores

class NeighborDistanceScoring(ScoringMethod):
    def __init__(self):
        super().__init__('Neighbor Distance Scoring')

    def calculate(self, snp, snps, reference, genes, tandem_repeats):
        sorted_snps = sorted(snps)
        snp_index = sorted_snps.index(snp)
        distance_previous, distance_next = None, None
        if snp_index > 0:
            distance_previous = snp - sorted_snps[snp_index - 1]
        if snp_index < len(sorted_snps) - 1:
            distance_next = sorted_snps[snp_index + 1] - snp
        min_distance = min(distance_previous, distance_next) if distance_previous is not None and distance_next is not None else distance_previous or distance_next
        return min_distance

class GCContentScoring(ScoringMethod):
    def __init__(self, window_size=100):
        super().__init__('GC Content Scoring')
        self.window_size = window_size

    def calculate(self, snp, snps, reference, genes, tandem_repeats):
        start = max(0, snp - self.window_size // 2)
        end = min(len(reference), snp + self.window_size // 2)
        sequence = reference[start:end]
        gc_count = sequence.count('G') + sequence.count('C')
        gc_content = gc_count / len(sequence) * 100
        return gc_content


class TandemRepeatProximityScoring(ScoringMethod):
    def __init__(self, window_size=100):
        super().__init__('Tandem Repeat Proximity Scoring')
        self.window_size = window_size

    def calculate(self, snp, snps, reference, genes, tandem_repeats):
        min_distance = float('inf')

        # Consider SNPs within the window size
        start = max(0, snp - self.window_size // 2)
        end = min(len(reference), snp + self.window_size // 2)

        for interval in tandem_repeats:
            repeat_start, repeat_end = interval
            if repeat_end < start:
                continue
            if repeat_start > end:
                break
            distance = min(abs(snp - repeat_start), abs(snp - repeat_end))
            min_distance = min(min_distance, distance)

        return min_distance if min_distance != float('inf') else None


class GeneProximityScoring(ScoringMethod):
    def __init__(self, window_size=100):
        super().__init__('Gene Proximity Scoring')
        self.window_size = window_size

    def calculate(self, snp, snps, reference, genes, tandem_repeats):
        min_distance = float('inf')

        # Consider SNPs within the window size
        start = max(0, snp - self.window_size // 2)
        end = min(len(reference), snp + self.window_size // 2)

        for gene in genes:
            gene_start, gene_end, gene_type = gene
            if gene_end < start:
                continue
            if gene_start > end:
                break

            distance = min(abs(snp - gene_start), abs(snp - gene_end))
            min_distance = min(min_distance, distance)

        return min_distance if min_distance != float('inf') else None


class MutationTypeScoring(ScoringMethod):
    def __init__(self):
        super().__init__('Mutation Type Scoring')

    def calculate(self, snp, snps, reference, genes, tandem_repeats):
        # Assumed that each gene info contains: [start_position, end_position, gene_type]
        # And gene_type can be 'synonymous', 'nonsynonymous', or 'other'

        mutation_type = 'other'

        for gene in genes:
            start, end, gene_type = gene
            if start <= snp <= end:
                mutation_type = gene_type
                break

        # Scoring based on the mutation type
        if mutation_type == 'synonymous':
            return 1
        elif mutation_type == 'nonsynonymous':
            return 2
        else:
            return 0
