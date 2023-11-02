from cansnpscorer.helpers.data_utils import neighbour_distance, gc_content, search, gene_checker


def NeighborDistanceScoring(snp_idx, snps, reference, correction, chromosome_length):
    """
    Wrapper around neighbour_distance for easier usage in ScoringManager.
    """
    current_snp = snps[snp_idx]

    # Logic for getting left and right neighbors from original code
    if snp_idx > 0 and snp_idx < len(snps) - 1:
        left = snps[snp_idx - 1][0]
        right = snps[snp_idx + 1][0]
    elif snp_idx == 0:
        left = 0
        right = snps[snp_idx + 1][0]
    elif snp_idx == len(snps) - 1:
        left = snps[snp_idx - 1][0]
        right = chromosome_length  # length of the chromosome

    # Apply correction for chromosome switching
    if current_snp[0] - correction > chromosome_length:
        correction = reference[chromosome_length]

    snp_pos = current_snp[0] - correction

    return neighbour_distance(snp_pos, left, right)


def GCContentScoring(snp_idx, snps, reference, correction, gc_goal, chromosome_length):
    """
    Wrapper around gc_content for easier usage in ScoringManager.
    """
    current_snp = snps[snp_idx]

    # Apply correction for chromosome switching
    if current_snp[0] - correction > chromosome_length:
        correction = reference[chromosome_length]

    snp_pos = current_snp[0] - correction

    return gc_content(snp_pos, reference[0], gc_goal)


def TandemRepeatProximityScoring(snp_idx, snps, trf, correction, chromosome, tandem_limit, chromosome_length):
    """
    Wrapper around search for tandem repeats for easier usage in ScoringManager.
    """
    current_snp = snps[snp_idx]

    # Apply correction for chromosome switching
    if current_snp[0] - correction > chromosome_length:
        correction = trf[chromosome]
        chromosome += 1

    snp_pos = current_snp[0] - correction

    return search(snp_pos, trf[chromosome], tandem_limit)[0]

def GeneProximityScoring(snp_idx, snps, ref_gb, correction, chromosome, gene_limit, chromosome_length):
    """
    Wrapper around search for gene proximity for easier usage in ScoringManager.
    """
    current_snp = snps[snp_idx]

    # Apply correction for chromosome switching
    if current_snp[0] - correction > chromosome_length:
        correction = ref_gb[chromosome]
        chromosome += 1

    snp_pos = current_snp[0] - correction

    gene_proximity = search(snp_pos, ref_gb[chromosome], gene_limit)

    # Adjust the gene score for aggregation
    adjusted_gene_score = (1 - gene_proximity[0])

    return adjusted_gene_score, gene_proximity



def MutationTypeScoring(snp_idx, snps, reference, gene_score, correction):
    """
    Wrapper around gene_checker for mutation type scoring for easier usage in ScoringManager.
    """
    current_snp = snps[snp_idx]

    # If gene_score[0] is 1, it means the SNP is not within or near any gene; thus mutation type is irrelevant
    if gene_score[0] != 0:
        gene_proximity = gene_score[1]
        gene_interval = gene_proximity[1]
        mutation_type_score = gene_checker(current_snp, reference, gene_interval, correction)
        return mutation_type_score
    else:
        return 0




def ScoringManager(i, snps, reference, ref_gb, trf, chromosome_length, correction, chromosome,
                   tandem_limit=50, gene_limit=0, gc_goal=0.5, modifier=[]):
    """
    Aggregates all the scoring functions and calculates the final score.
    """
    # Individual Scoring Functions
    neighbour_score = NeighborDistanceScoring(i, snps, reference, correction, chromosome_length)
    gc_score = GCContentScoring(i, snps, reference, correction, gc_goal, chromosome_length)
    tandem_score = TandemRepeatProximityScoring(i, snps, trf, correction, chromosome, tandem_limit, chromosome_length)
    gene_score = GeneProximityScoring(i, snps, ref_gb, correction, chromosome, gene_limit, chromosome_length)
    mutation_score = MutationTypeScoring(i, snps, reference, gene_score, correction)

    # Mapping scores to their names
    scores_dict = {
        'neighbour_score': neighbour_score,
        'gc_score': gc_score,
        'tandem_score': tandem_score,
        'gene_score': gene_score[0],
        'mutation_score': mutation_score
    }

    # Apply the modifier factor to selected scores and multiply all by 100
    for score_name, value in scores_dict.items():
        factor = 2 if score_name in modifier else 1
        scores_dict[score_name] = round(value * factor * 100)

    # Aggregate the scores
    score = list(scores_dict.values())

    total = sum(score)

    # The final score array consists of SNP position, SNP variant, total score, and individual scores.
    final_score = [snps[i][0], snps[i][1], total] + score

    return final_score
