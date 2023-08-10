import sys


class Parameters:

    # Set default values for parameters
    run_number = 0
    num_of_runs = 1  # should be 1000
    sample_size = 10  # should be 1000
    n_estimators = 100  # should be ??
    sequence_length = 1000  # should be 1000
    snp_prob = 0.01  # should be ??

    inversion_prob = 0  # should be 0.3
    inversion_length_min = 5  # should be 5
    inversion_length_max = 5  # should be 30

    ir_force_prob = 1  # should be ??
    ir_inversion_prob = 0.3  # should be 0.3; If ir exists what is the probability it'll be inverted
    ir_arm_length = 10  # should be ??
    ir_spacer_length = 5  # should be ??

    def __init__(self):
        # Read command-line arguments and update parameters if provided
        # reminder that sys.argv[0] is the script name
        if len(sys.argv) == 1:
            return  # use all default parameters
        if len(sys.argv) < 13:
            raise ValueError("Not enough input parameters. Expecting {} got {}.".format(13, len(sys.argv)))

        self.run_number = int(Parameters.parameter_clean_up(sys.argv[1]))
        self.num_of_runs = int(Parameters.parameter_clean_up(sys.argv[2]))
        self.sample_size = int(Parameters.parameter_clean_up(sys.argv[3]))
        self.n_estimators = int(Parameters.parameter_clean_up(sys.argv[4]))
        self.sequence_length = int(Parameters.parameter_clean_up(sys.argv[5]))
        self.snp_prob = int(Parameters.parameter_clean_up(sys.argv[6]))

        self.inversion_prob = float(Parameters.parameter_clean_up(sys.argv[7]))
        self.inversion_length_min = int(Parameters.parameter_clean_up(sys.argv[8]))
        self.inversion_length_max = int(Parameters.parameter_clean_up(sys.argv[9]))

        self.ir_force_prob = float(Parameters.parameter_clean_up(sys.argv[10]))
        self.ir_inversion_prob = float(Parameters.parameter_clean_up(sys.argv[11]))
        self.ir_arm_length = int(Parameters.parameter_clean_up(sys.argv[12]))
        self.ir_spacer_length = int(Parameters.parameter_clean_up(sys.argv[13]))

        if self.inversion_prob > 0 and (self.ir_force_prob > 0 or self.ir_inversion_prob > 0):
            raise ValueError("Only normal inversion OR IR inversion is supported in a single run")

    @staticmethod
    def parameter_clean_up(parameter_index):
        return sys.argv[parameter_index].split('#', 1)[0].strip()