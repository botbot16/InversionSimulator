import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from util import Util


class RandomForest:

    def __init__(self, n_estimators):
        self.n_estimators = n_estimators

    @staticmethod
    def generate_features(original, mutated):
        sample_size = len(original)
        if sample_size != len(mutated):
            raise AssertionError("Array sizes should be equal.")

        sequence_length = []
        match_score = []
        is_inverted_repeat = []
        for i in range(0, sample_size):
            sequence_length.append(len(original))
            match_score.append(Util.calc_match_score_pair(original[i], mutated[i]))
            is_inverted_repeat.append(Util.has_inverted_repeat(original[i]))

        return np.column_stack((sequence_length, match_score, is_inverted_repeat))

    @staticmethod
    def generate_features_leaves(leaf1, leaf2, leaf3):
        if not len(leaf1) == len(leaf2) == len(leaf3):
            raise AssertionError("Array sizes should be equal")
        sample_size = len(leaf1)
        sequence_length = []
        match_score = []
        is_inverted_repeat = []
        for i in range(0, sample_size):
            if not len(leaf1[i]) == len(leaf2[i]) == len(leaf3[i]):
                raise AssertionError("Leaf sizes should be equal")
            sequence_length.append(len(leaf1))
            match_score.append(Util.calc_match_score_pair(leaf1[i], leaf2[i]) +
                               Util.calc_match_score_pair(leaf1[i], leaf3[i]) +
                               Util.calc_match_score_pair(leaf2[i], leaf3[i]))
            is_inverted_repeat.append(Util.has_inverted_repeat_triple(leaf1[i], leaf2[i], leaf3[i]))

        return np.column_stack((sequence_length, match_score, is_inverted_repeat))

    def learn(self, data, target):
        # split into train and test sets:
        data_train, data_test, target_train, target_test = train_test_split(data, target, test_size=0.3)
        rf_classifier = RandomForestClassifier(n_estimators=self.n_estimators)
        rf_classifier.fit(data_train, target_train)
        model_predictions = rf_classifier.predict(data_test)
        # print(model_predictions)
        score = accuracy_score(model_predictions, target_test)
        return score



