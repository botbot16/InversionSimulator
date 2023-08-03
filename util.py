import sys
import pandas
from Bio import pairwise2

NUCLEOTIDES = ('A', 'C', 'G', 'T')
NUCLEOTIDES_INVERSION = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


class Util:

    @staticmethod
    def read_parameters():
        # Set default values for parameters
        sequence_length = 1000  # should be 1000
        inversion_prob = 0.3  # should be 0.3
        inversion_length_min = 5  # should be 5
        inversion_length_max = 30  # should be 30
        ir_force_prob = 0.3  # should be ??
        ir_inversion_prob = 0.3  # should be 0.3; If ir exists what is the probability it'll be inverted
        ir_arm_length = 7  # should be ??
        ir_spacer_length = 4  # should be ??
        snp_prob = 0.001  # should be ??
        n_estimators = 100  # should be ??
        sample_size = 1000  # should be 1000

        # Read command-line arguments and update parameters if provided
        # reminder that sys.argv[0] is the script name
        if len(sys.argv) >= 2:
            sequence_length = int(sys.argv[1])
        if len(sys.argv) >= 5:
            inversion_prob = float(sys.argv[2])
            inversion_length_min = int(sys.argv[3])
            inversion_length_max = int(sys.argv[4])
        if len(sys.argv) >= 9:
            ir_force_prob = float(sys.argv[5])
            ir_inversion_prob = float(sys.argv[6])
            ir_arm_length = int(sys.argv[7])
            ir_spacer_length = int(sys.argv[8])
        if len(sys.argv) >= 10:
            snp_prob = int(sys.argv[9])
        if len(sys.argv) >= 11:
            n_estimators = int(sys.argv[10])
            sample_size = int(sys.argv[11])

        # Return the parameters
        return sequence_length, inversion_prob, inversion_length_min, inversion_length_max, ir_force_prob, \
               ir_inversion_prob, ir_arm_length, ir_spacer_length, snp_prob, n_estimators, sample_size

    @staticmethod
    def inversion(sequence):
        return ''.join(NUCLEOTIDES_INVERSION[nuc] for nuc in sequence[::-1])

    @staticmethod
    def force_inverted_repeat(sequence, ir_arm_length, ir_spacer_length):
        if len(sequence) != 2 * ir_arm_length + ir_spacer_length:
            raise AssertionError("sequence length doesn't match")
        sequence_with_ir = sequence[:-ir_arm_length] + Util.inversion(sequence[:ir_arm_length])
        # print("force_inverted_repeat with " + str(ir_arm_length) + " and " + str(ir_spacer_length))
        # print(sequence)
        # print(sequence_with_ir)
        return sequence_with_ir

    @staticmethod
    def calc_match_score_pair(seq1, seq2):
        alignments = pairwise2.align.globalxx(seq1, seq2)
        best_alignment = alignments[0]  # Get the best alignment
        alignment_score = best_alignment[2]
        return alignment_score

    @staticmethod
    def has_inverted_repeat(sequence):
        inverted_repeats = {}
        length = len(sequence)
        for ir_size in range(7, 8):
            for i in range(length):
                subsequence = sequence[i:i + ir_size]
                if subsequence in inverted_repeats and inverted_repeats[subsequence] + ir_size + 1 < i \
                        and inverted_repeats[subsequence] + ir_size + 4 == i:
                    # print("IR found")
                    # print(sequence)
                    # print(subsequence + " at pos " + str(i))
                    # print(Util.inversion(subsequence) + " at pos " + str(inverted_repeats[subsequence]))
                    return 1
                inverted_repeats[Util.inversion(subsequence)] = i

        # print("IR not found")
        return 0  # no reverse compliments :(

    @staticmethod
    def has_inverted_repeat_triple(seq1, seq2, seq3):
        return Util.has_inverted_repeat(seq1) or Util.has_inverted_repeat(seq2) or Util.has_inverted_repeat(seq3)
