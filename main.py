import os
from parameters import Parameters
from sequence_wrapper import SequenceWrapper
from random_forest import RandomForest
import datetime
import numpy as np


def get_num_of_leaves_as_string():
    if num_of_leaves == 2:
        return "two"
    elif num_of_leaves == 3:
        return "three"
    elif num_of_leaves == 4:
        return "four"


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


def generate_leaves(parameters, num_of_leaves):
    leaf1 = []
    leaf2 = []
    leaf3 = []
    leaf4 = []
    target = []
    for i in range(0, parameters.sample_size):
        sequence_wrapper = SequenceWrapper(parameters.sequence_length, parameters.ir_force_prob,
                                           parameters.ir_arm_length, parameters.ir_spacer_length)
        sequence_wrapper.mutate_n(num_of_leaves, parameters.ir_inversion_prob, parameters.ir_arm_length,
                                  parameters.ir_spacer_length, parameters.inversion_prob,
                                  parameters.inversion_length_min, parameters.inversion_length_max,
                                  parameters.snp_prob)
        leaf1.append(sequence_wrapper.leaf1)
        leaf2.append(sequence_wrapper.leaf2)
        leaf3.append(sequence_wrapper.leaf3)
        leaf4.append(sequence_wrapper.leaf4)
        target.append(sequence_wrapper.was_inverted)

    random_forest = RandomForest(parameters.n_estimators)
    data = random_forest.generate_features_leaves(parameters, num_of_leaves, leaf1, leaf2, leaf3, leaf4)
    return random_forest.learn(data, target)


def print_results(file_path, scores):
    os.makedirs(file_path, exist_ok=True)
    file_path += run_number + ".txt"
    with open(file_path, 'w') as file:
        file.write("median: {} min: {}, max: {}\nfull scores: {}".format(
            np.median(scores), min(scores), max(scores),
            scores))


if __name__ == "__main__":
    start_time = datetime.datetime.now()
    print("Run start time: ", start_time)
    parameters = Parameters()
    root_and_leaf_scores = []
    leaf_scores = []
    formatted_date = start_time.strftime("%Y%m%d_%H%M")

    for i in range(parameters.num_of_runs):
        root_and_leaf_scores.append(root_and_leaf(parameters))

    path_prefix = formatted_date + "/"
    path_prefix += "IR/" if parameters.inversion_prob == 0 else "inversion/"
    run_number = str(parameters.run_number)

    print_results(path_prefix + "root_and_leaf/", root_and_leaf_scores)

    num_of_leaves = parameters.num_of_leaves
    string_num_of_leaves = get_num_of_leaves_as_string()

    for i in range(parameters.num_of_runs):
        leaf_scores.append(generate_leaves(parameters, num_of_leaves))

    print_results(path_prefix + string_num_of_leaves + "_leaves/", leaf_scores)

    end_time = datetime.datetime.now()
    print("Run end time:", end_time)
    print("Run length is:", end_time - start_time)
