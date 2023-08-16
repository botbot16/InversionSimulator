import random
import util
from util import Util

class SequenceWrapper:

    def __init__(self, sequence_length, ir_prob, ir_arm_length, ir_spacer_length):
        self.force_ir_position = None
        self.mutated_sequence, self.leaf1, self.leaf2, self.leaf3, self.leaf4 = None, None, None, None, None
        self.was_inverted = False

        # sequence_length random nucleotides with equal probability
        generated_sequence = ''.join(random.choice(util.NUCLEOTIDES) for _ in range(sequence_length))
        inverted_repeat_length = ir_spacer_length + 2 * ir_arm_length
        if random.random() >= ir_prob or inverted_repeat_length > sequence_length:  # no inverted repeat
            self.original_sequence = generated_sequence
        else:
            force_ir_position = random.randint(0, sequence_length - inverted_repeat_length)
            subsequence_with_ir = Util.force_inverted_repeat(
                generated_sequence[force_ir_position:force_ir_position + inverted_repeat_length],
                ir_arm_length, ir_spacer_length)
            result_sequence = generated_sequence[:force_ir_position] + subsequence_with_ir + \
                              generated_sequence[force_ir_position + inverted_repeat_length:]
            self.original_sequence = result_sequence
            self.force_ir_position = force_ir_position

    def mutate(self, ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob, inversion_length_min,
               inversion_length_max, snp_prob):
        sequence_length = len(self.original_sequence)
        result_sequence = self.original_sequence
        inversion_length = random.randint(inversion_length_min, inversion_length_max)
        if sequence_length < inversion_length:
            raise ValueError("sequence length shorter than requested inversion.")

        was_inverted = False
        # adding an ir_inversion
        if self.force_ir_position is not None and random.random() < ir_inversion_prob:
            was_inverted = True
            inversion_start_position = self.force_ir_position + ir_arm_length
            inversion_end_position = inversion_start_position + ir_spacer_length
            inverted_sequence = Util.inversion(
                result_sequence[inversion_start_position:inversion_end_position])
            result_sequence = result_sequence[:inversion_start_position] + inverted_sequence + \
                              result_sequence[inversion_end_position:]

        # adding an inversion
        if random.random() < inversion_prob:
            was_inverted = True
            position = random.randint(0, sequence_length - inversion_length)
            inverted_sequence = Util.inversion(result_sequence[position:position + inversion_length])
            result_sequence = result_sequence[:position] + inverted_sequence + result_sequence[
                                                                               position + inversion_length:]

        # SNPing every nucleotide with prob snp_prob (with a 1/4 chance for each possible nucleotide)
        if snp_prob > 0:
            tmp_sequence = ""
            for n in result_sequence:
                if random.random() < snp_prob:
                    tmp_sequence += random.choice(util.NUCLEOTIDES)
                else:
                    tmp_sequence += n
            result_sequence = tmp_sequence

        return result_sequence, was_inverted

    def mutate_one(self, ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob, inversion_length_min,
                   inversion_length_max, snp_prob):
        self.mutated_sequence, self.was_inverted = self.mutate(ir_inversion_prob, ir_arm_length, ir_spacer_length,
                                                               inversion_prob, inversion_length_min,
                                                               inversion_length_max, snp_prob)

    def mutate_n(self, n, ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob, inversion_length_min,
                 inversion_length_max, snp_prob):
        self.leaf1, was_inverted1 = self.mutate(ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob,
                                                inversion_length_min, inversion_length_max, snp_prob)
        self.leaf2, was_inverted2 = self.mutate(ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob,
                                                inversion_length_min, inversion_length_max, snp_prob)
        was_inverted = was_inverted1 or was_inverted2
        if n > 2:
            self.leaf3, was_inverted3 = self.mutate(ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob,
                                                    inversion_length_min, inversion_length_max, snp_prob)
            was_inverted = was_inverted or was_inverted3
        if n > 3:
            self.leaf4, was_inverted4 = self.mutate(ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob,
                                                    inversion_length_min, inversion_length_max, snp_prob)
            was_inverted = was_inverted or was_inverted4
        self.was_inverted = was_inverted

