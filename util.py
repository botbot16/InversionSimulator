from Bio import pairwise2

NUCLEOTIDES = ('A', 'C', 'G', 'T')
NUCLEOTIDES_INVERSION = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


class Util:

    @staticmethod
    def inversion(sequence):
        return ''.join(NUCLEOTIDES_INVERSION[nuc] for nuc in sequence[::-1])

    @staticmethod
    def force_inverted_repeat(sequence, ir_arm_length, ir_spacer_length):
        if len(sequence) != 2 * ir_arm_length + ir_spacer_length:
            raise ValueError("sequence length doesn't match")
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
    def has_inverted_repeat(parameters, sequence):
        inverted_repeats = {}
        length = len(sequence)
        for i in range(length):
            subsequence = sequence[i:i + parameters.ir_arm_length]
            # and inverted_repeats[subsequence] + Util.ir_arm_length + 1 < i \
            if subsequence in inverted_repeats \
                    and inverted_repeats[subsequence] + parameters.ir_arm_length + parameters.ir_spacer_length == i:
                # print(sequence)
                # print(subsequence + " at pos " + str(i))
                # print(Util.inversion(subsequence) + " at pos " + str(inverted_repeats[subsequence]))
                return 1
            inverted_repeats[Util.inversion(subsequence)] = i

        # print("IR not found")
        return 0  # no reverse compliments :(

    @staticmethod
    def has_inverted_repeat_triple(parameters, seq1, seq2, seq3):
        return Util.has_inverted_repeat(parameters, seq1) or Util.has_inverted_repeat(parameters, seq2) or \
               Util.has_inverted_repeat(parameters, seq3)
