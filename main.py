import os
from parameters import Parameters
from sequence_wrapper import SequenceWrapper
from random_forest import RandomForest
import datetime
import numpy as np


def root_and_leaf(parameters):
    original = []
    mutated = []
    target = []
    for i in range(0, parameters.sample_size):
        sequence_wrapper = SequenceWrapper(parameters.sequence_length, parameters.ir_force_prob,
                                           parameters.ir_arm_length, parameters.ir_spacer_length)
        sequence_wrapper.mutate_one(parameters.ir_inversion_prob, parameters.ir_arm_length, parameters.ir_spacer_length,
                                    parameters.inversion_prob, parameters.inversion_length_min,
                                    parameters.inversion_length_max, parameters.snp_prob)
        original.append(sequence_wrapper.original_sequence)
        mutated.append(sequence_wrapper.mutated_sequence)
        target.append(sequence_wrapper.was_inverted)

    random_forest = RandomForest(parameters.n_estimators)
    data = random_forest.generate_features(parameters, original, mutated)
    return random_forest.learn(data, target)

def three_leaves(parameters):
    leaf1 = []
    leaf2 = []
    leaf3 = []
    target = []
    for i in range(0, parameters.sample_size):
        sequence_wrapper = SequenceWrapper(parameters.sequence_length, parameters.ir_force_prob,
                                           parameters.ir_arm_length, parameters.ir_spacer_length)
        sequence_wrapper.mutate_three(parameters.ir_inversion_prob, parameters.ir_arm_length,
                                      parameters.ir_spacer_length, parameters.inversion_prob,
                                      parameters.inversion_length_min, parameters.inversion_length_max,
                                      parameters.snp_prob)
        leaf1.append(sequence_wrapper.leaf1)
        leaf2.append(sequence_wrapper.leaf2)
        leaf3.append(sequence_wrapper.leaf3)
        target.append(sequence_wrapper.was_inverted)

    random_forest = RandomForest(parameters.n_estimators)
    data = random_forest.generate_features_leaves(parameters, leaf1, leaf2, leaf3)
    return random_forest.learn(data, target)


if __name__ == "__main__":
    parameters = Parameters()
    root_and_leaf_scores = []
    three_leaves_scores = []
    formatted_date = datetime.datetime.now().strftime("%Y%m%d_%H%M")

    for i in range(parameters.num_of_runs):
        root_and_leaf_scores.append(root_and_leaf(parameters))

    path_prefix = formatted_date + "/"
    path_prefix += "IR/" if parameters.inversion_prob == 0 else "inversion/"
    run_number = str(parameters.run_number)

    file_path = path_prefix + "root_and_leaf/"
    os.makedirs(file_path)
    file_path += run_number + ".txt"
    with open(file_path, 'w') as file:
        file.write("root_and_leaf statistics:\nmedian: {} min: {}, max: {}\nfull scores: {}".format(
                   np.median(root_and_leaf_scores), min(root_and_leaf_scores), max(root_and_leaf_scores),
                   root_and_leaf_scores))

    for i in range(parameters.num_of_runs):
        three_leaves_scores.append(three_leaves(parameters))

    file_path = path_prefix + "three_leaves/"
    os.makedirs(file_path)
    file_path += run_number + ".txt"
    with open(file_path, 'w') as file:
        file.write("median statistics:\nmedian: {} min: {}, max: {}\nfull scores: {}".format(
                   np.median(three_leaves_scores), min(three_leaves_scores), max(three_leaves_scores),
                   three_leaves_scores))
