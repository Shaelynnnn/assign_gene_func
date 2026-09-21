import numpy as np
from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman-Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    # get the n and m
    len1 = len(seq1)
    len2 = len(seq2)
    # creat a 2D matric (score matric)
    score_matrix = np.zeros((len1 + 1, len2 + 1))
    # initialise the row & col;
    # row;
    i = 1
    while i <= len1:
        score_matrix[i, 0] = score_matrix[i - 1, 0] + scoring_function(seq1[i - 1], '-')
        i += 1
    # col;
    j = 1
    while j <= len2:
        score_matrix[0, j] = score_matrix[0, j - 1] + scoring_function('-', seq2[j - 1])
        j += 1
    
    # use dp to complete the remain cells;
    i = 1
    while i <= len1:
        j = 1
        while j <= len2:
            align = score_matrix[i - 1, j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            insertion = score_matrix[i - 1, j] + scoring_function(seq1[i - 1], "-")
            deletion = score_matrix[i, j - 1] + scoring_function("-", seq2[j - 1])
            # get the current max;
            score_matrix[i, j] = max(align, deletion, insertion)

            j += 1
        i += 1
        
    # traceback;
    # get the final align mark
    score = score_matrix[len1, len2]
    # start point;
    i = len1
    j = len2
    # prepare the results
    seq1_align = []
    seq2_align = []
    
    # check current score;
    while i != 0 or j != 0:
        align = -np.inf
        insertion = -np.inf
        deletion = -np.inf
        
        if i > 0 and j > 0:
            align = score_matrix[i - 1, j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
        if i > 0:
            insertion = score_matrix[i - 1, j] + scoring_function(seq1[i - 1], "-")
        if j > 0:
            deletion = score_matrix[i, j - 1] + scoring_function("-", seq2[j - 1])
        
        # dignose
        if score_matrix[i, j] == align:
            seq1_align.append(seq1[i - 1])
            seq2_align.append(seq2[j - 1])
            i -= 1
            j -= 1
        # up
        elif score_matrix[i, j] == insertion:
            seq1_align.append(seq1[i - 1])
            seq2_align.append("-")
            i -= 1
        # left
        elif score_matrix[i, j] == deletion:
            seq1_align.append("-")
            seq2_align.append(seq2[j - 1])
            j -= 1
    
    seq1_align.reverse()
    seq2_align.reverse()
    
    seq1_result = "".join(seq1_align)
    seq2_result = "".join(seq2_align)
    
    return (seq1_result, seq2_result, float(score))


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

# use the BLOSUM62
def scoring_function_blosum62(aa_i, aa_j):
    if aa_i == '-' or aa_j == '-':
        return -5
    
    return float(blosum62[aa_i, aa_j])
