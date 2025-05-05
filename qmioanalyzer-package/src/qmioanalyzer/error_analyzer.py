import numpy as np

class ErrorAnalyzer:
    '''
    Class to compute the number of errors given an array of bits.

    To compute the number of errors we consider 2 cases:

    1. Repetition_period > 500 mus. 

        In this case the system is expected to go back to |0> state
        in each iteration, thus the definition of error will only depend 
        on the initialization of the state in the circuit. 
        - If an X gate is applied: 
            The state initializes to |1> -> measuring 0 is an error.
        - If nothing is applied:
            The state initializes to |0> -> measuring 1 is an error.
    
    2. Repetition period < 500 mus.

        In this case the system doesn't recover after applying an X gate.
        The expected circuit would consist in applying an X gate, waiting 
        a certain delay time and measuring the qubits, executed for each 
        *repetition_period* s.

        The definition of the errors here would be given by the previous 
        measurement. As an example:

        A measurement in the iteration i gives 1101, for i+1 we expect 0010;
        instead 0000 is measured, registering 1 error.

        A classification between the possible errors is made here:

        A. True errors: A 0 measurement is expected, instead we get a 1.
        B. False errors: A 1 measurement is expected, instead we get a 0.
        
    '''

    def __init__(self, repetition_period, init_state, measurements, verbose=0):
        """
        Parameters:
        - repetition_period: float, time between measurements in s.
        - init_state: initialization of the state : 
                      0 -> if nothing is applied
                      1 -> if an x gate is applied
        - measurements: array of measured eigenvalues (of type string).
        - verbose : parameter to control the explanation given of the analysis.
                     0 -> no explanation.
                     1 -> short  summary of the returned dictionary
                     2 -> extended summary of the result
        """

        self.repetition_period = repetition_period
        self.init_state = init_state
        self.measurements = measurements
        self.errors = {} 
        self.analysis_type = None  # Store the type of analysis used

        self.verbose = verbose  # verbosity level
    
    def analyze_errors(self):
        '''
        Method to compute the number of errors according to each case.
        '''
        if self.repetition_period >= 0.0005:
            self.analysis_type = "high_period"
            self._analyze_high_period()
        else:
            self.analysis_type = "low_period"
            self._analyze_low_period()
    
    def _analyze_high_period(self):
        '''
        Method to compute the number of errors when there is enough time
        between measurements for the state to go back to |0>, i.e.,
        when repetition_period >= 500 mus. 
        '''
        
        expected_value = '0' if self.init_state == 0 else '1'

        num_qubits = len(self.measurements[0])
        errors_per_qubit = np.zeros(num_qubits)


        # we initialize the arrays to save the results
        self.errors["high_period_errors"] = {
            "error_counts": [],
            "times": [],  # array to save the corresponding time to each measurement
            "errors_per_qubit": errors_per_qubit
        }

        for i, measurement in enumerate(self.measurements):
            error_count = 0  
            for qubit_index, bit in enumerate(measurement):
                if bit != expected_value:
                    errors_per_qubit[qubit_index] += 1  # Contamos el error para ese qubit
                    error_count += 1  


            measurement_time = i * self.repetition_period # we compute the time from the index and the repetition_period
            self.errors["high_period_errors"]["error_counts"].append(error_count)
            self.errors["high_period_errors"]["times"].append(measurement_time)
        
        self.errors["high_period_errors"]["states"] = [state for state in self.measurements]
        self.errors["high_period_errors"]["errors_per_qubit"] = errors_per_qubit.tolist()
    
    def _analyze_low_period(self):
        '''
        Method to compute the number of errors when X gates are applied 
        continuously, i.e., when repetition_period < 500 mus.

        The number of errors is defined from the previous measurement.

        A classification of the errors is made here in True and False errors,
        according to the already given definitions.
        '''

        expected_measurements = self._compute_expected_sequences()
        num_qubits = len(self.measurements[0])
        errors_per_qubit = np.zeros(num_qubits)

        # we initialize the arrays to save the results 
        self.errors["low_period_errors"] = {
                    "true_errors": [],
                    "false_errors": [],
                    "total_errors": [],
                    "times": [],
                    "errors_per_qubit": errors_per_qubit
                }


        for i, (expected, measured) in enumerate(zip(expected_measurements, self.measurements[1:])):
            
            true_errors = 0
            false_errors = 0

            for qubit_index, (e, m) in enumerate(zip(expected, measured)):
                #print(e,m)
                if e == '0' and m == '1':
                    true_errors += 1
                    errors_per_qubit[qubit_index] += 1
                elif e == '1' and m == '0':
                    false_errors += 1
                    errors_per_qubit[qubit_index] += 1

            measurement_time = i * self.repetition_period
            
            self.errors["low_period_errors"]["true_errors"].append(true_errors)
            self.errors["low_period_errors"]["false_errors"].append(false_errors)
            self.errors["low_period_errors"]["total_errors"].append(true_errors + false_errors)
            self.errors["low_period_errors"]["times"].append(measurement_time)

        self.errors["low_period_errors"]["states"] = self.measurements[1:]
        self.errors["low_period_errors"]["expected"] = expected_measurements
        self.errors["low_period_errors"]["errors_per_qubit"] = errors_per_qubit.tolist()


    
    def _compute_expected_sequences(self):
        '''
        Method to compute the expected values of the bits from 
        the previous ones simulating the application of the X gate.
        '''
        expected_sequences = []
        
        previous = self.measurements[0]
        for i in range(1, len(self.measurements)):  
            expected = ''.join('1' if b == '0' else '0' for b in previous)  # simulating an X gate
            expected_sequences.append(expected)
            previous = self.measurements[i]  # we use the measurement as a reference 
        return expected_sequences
    
    def get_errors(self):
        '''
        Returns the computed errors as a dictionary.
        '''
        self.analyze_errors()

        if self.verbose == 1:
            self._get_short_summary()
        elif self.verbose == 2:
            self._get_extended_summary()

        return self.errors

    def _get_extended_summary(self):
        '''
        Returns an extended summary of the error analysis performed.
        '''
        if self.analysis_type is None:
            print("No error analysis has been performed yet.")
        
        if self.analysis_type == "high_period":
            print("The analysis was performed using the *high repetition period* approach ")
            print("(repetition_period ≥ 500 μs). The error dictionary contains:\n")
            print(f"- Key: 'high_period_errors' → Returns two arrays with {len(self.errors['high_period_errors']['error_counts'])} measurements:")
            print("  - 'error_counts': Array with the number of errors associated with each measurement.")
            print("  - 'times': Array with the times of each measurement (calculated as i * repetition_period).")
            print("\n")
        else:
            print("The analysis was performed using the *low repetition period* approach ")
            print("(repetition_period < 500 μs). The error dictionary contains:\n")
            print(f"- Key: 'low_period_errors' → A dictionary with {len(self.errors['low_period_errors']['true_errors'])} measurements with the following subkeys:")
            print("  - 'true_errors': Array with the number of true errors (0 expected, 1 measured) for each measurement.")
            print("  - 'false_errors': Array with the number of false errors (1 expected, 0 measured) for each measurement.")
            print("  - 'total_errors': Array with the total number of errors for each measurement.")
            print("  - 'times': Array with the times of each measurement (calculated as i * repetition_period).")
            print("\n")
        
    def _get_short_summary(self):
        if self.analysis_type == "high_period":
            print("Analysis mode: High Period")
            print("Stored in dictionary under key: 'high_period_errors' with two arrays: ")
            print(" - 'error_counts': Number of errors for each measurement.")
            print(" - 'times': Time of each measurement.")
            print("\n")
        else:
            print("Analysis mode: Low Period")
            print(f"Stored in dictionary under key: 'low_period_errors' with the following arrays:")
            print(" - 'true_errors': Number of true errors (0 expected, 1 measured).")
            print(" - 'false_errors': Number of false errors (1 expected, 0 measured).")
            print(" - 'total_errors': Total number of errors.")
            print(" - 'times': Time of each measurement.")
            print("\n")
        



