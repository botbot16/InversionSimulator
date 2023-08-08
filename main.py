from util import Util
from sequence_wrapper import SequenceWrapper
from random_forest import RandomForest
import numpy as np


def root_and_leaf():
    sequence_length, inversion_prob, inversion_length_min, inversion_length_max, ir_force_prob, ir_inversion_prob, \
    ir_arm_length, ir_spacer_length, snp_prob, n_estimators, sample_size = Util.read_parameters()
    original = []
    mutated = []
    target = []
    for i in range(0, sample_size):
        sequence_wrapper = SequenceWrapper(sequence_length, ir_force_prob, ir_arm_length, ir_spacer_length)
        sequence_wrapper.mutate_one(ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob,
                                    inversion_length_min, inversion_length_max, snp_prob)
        original.append(sequence_wrapper.original_sequence)
        mutated.append(sequence_wrapper.mutated_sequence)
        target.append(sequence_wrapper.was_inverted)

    random_forest = RandomForest(n_estimators)
    data = random_forest.generate_features(original, mutated)
    return random_forest.learn(data, target)

def three_leaves():
    sequence_length, inversion_prob, inversion_length_min, inversion_length_max, ir_force_prob, ir_inversion_prob, \
    ir_arm_length, ir_spacer_length, snp_prob, n_estimators, sample_size = Util.read_parameters()
    leaf1 = []
    leaf2 = []
    leaf3 = []
    target = []
    for i in range(0, sample_size):
        sequence_wrapper = SequenceWrapper(sequence_length, ir_force_prob, ir_arm_length, ir_spacer_length)
        sequence_wrapper.mutate_three(ir_inversion_prob, ir_arm_length, ir_spacer_length, inversion_prob,
                                      inversion_length_min, inversion_length_max, snp_prob)
        leaf1.append(sequence_wrapper.leaf1)
        leaf2.append(sequence_wrapper.leaf2)
        leaf3.append(sequence_wrapper.leaf3)
        target.append(sequence_wrapper.was_inverted)

    random_forest = RandomForest(n_estimators)
    data = random_forest.generate_features_leaves(leaf1, leaf2, leaf3)
    return random_forest.learn(data, target)


if __name__ == "__main__":
    root_and_leaf_scores = []
    three_leaves_scores = []
    for i in range(1000):
        root_and_leaf_scores.append(root_and_leaf())
    file_path = "root_and_leaf.txt"
    with open(file_path, 'w') as file:
        file.write("root_and_leaf statistics:\nmedian: {} min: {}, max: {}\nfull scores: {}".format(
                   np.median(root_and_leaf_scores), min(root_and_leaf_scores), max(root_and_leaf_scores),
                   root_and_leaf_scores))
    for i in range(1000):
        three_leaves_scores.append(three_leaves())
    file_path = "three_leaves.txt"
    with open(file_path, 'w') as file:
        file.write("median statistics:\nmedian: {} min: {}, max: {}\nfull scores: {}".format(
                   np.median(three_leaves_scores), min(three_leaves_scores), max(three_leaves_scores),
                   three_leaves_scores))
