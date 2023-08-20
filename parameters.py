import sys


class Parameters:

    # Set default values for parameters
    run_number = 0
    num_of_runs = 1  # should be 1000
    sample_size = 10  # should be 1000
    n_estimators = 10  # should be 10
    num_of_leaves = 4  # should be 2, 3, 4
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
        # reminder that sys.argv[0] is the script name, sys.argv[1] is argument string
        if len(sys.argv) > 1:
            arguments = sys.argv[1].split('\n')
            if len(arguments) < 13:
                raise ValueError("Not enough input parameters in sys.argv[1]. Expecting {} got {}."
                                 .format(13, len(arguments)))

            self.run_number = int(Parameters.parameter_clean_up(arguments[0]))
            self.num_of_runs = int(Parameters.parameter_clean_up(arguments[1]))
            self.sample_size = int(Parameters.parameter_clean_up(arguments[2]))
            self.n_estimators = int(Parameters.parameter_clean_up(arguments[3]))
            self.num_of_leaves = int (Parameters.parameter_clean_up(argument4))
            self.sequence_length = int(Parameters.parameter_clean_up(arguments[5]))
            self.snp_prob = float(Parameters.parameter_clean_up(arguments[6]))

            self.inversion_prob = float(Parameters.parameter_clean_up(arguments[7]))
            self.inversion_length_min = int(Parameters.parameter_clean_up(arguments[8]))
            self.inversion_length_max = int(Parameters.parameter_clean_up(arguments[9]))

            self.ir_force_prob = float(Parameters.parameter_clean_up(arguments[10]))
            self.ir_inversion_prob = float(Parameters.parameter_clean_up(arguments[11]))
            self.ir_arm_length = int(Parameters.parameter_clean_up(arguments[12]))
            self.ir_spacer_length = int(Parameters.parameter_clean_up(arguments[13]))

        if self.inversion_prob > 0 and (self.ir_force_prob > 0 or self.ir_inversion_prob > 0):
            raise ValueError("Only normal inversion OR IR inversion is supported in a single run")

        print("Starting a run with the following parameters:")
        print(', '.join([str(self.run_number), str(self.num_of_runs), str(self.sample_size), str(self.n_estimators),
                        str(self.sequence_length), str(self.snp_prob), str(self.inversion_prob),
                        str(self.inversion_length_min), str(self.inversion_length_max), str(self.ir_force_prob),
                        str(self.ir_inversion_prob), str(self.ir_arm_length), str(self.ir_spacer_length)]))

    @staticmethod
    def parameter_clean_up(parameter):
        return parameter.split('#', 1)[0].strip()
