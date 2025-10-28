from copy import error
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations

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

    def __init__(self, repetition_period,  measurements, threshold = 0.00015, init_state = None, 
                 verbose = 0, center_time=None, window=None, external_error_probs=None):
        """
        Parameters:
        - repetition_period: float, time between measurements in s.
        - init_state: initialization of the state : 
                      0 -> if nothing is applied
                      1 -> if an x gate is applied
                      None -> if the repetition_period is less than 500 mus
        - measurements: array of measured eigenvalues (of type string).
        - verbose : parameter to control the explanation given of the analysis.
                     0 -> no explanation.
                     1 -> short  summary of the returned dictionary
                     2 -> extended summary of the result
        - center_time: center time for windowing when a long measurement is done.
        - window: time window for selection
        - external_error_probs: array-like, optional. External probabilities for each qubit.
                           If provided, these will be used instead of computed probabilities
                           for theoretical histogram calculation.

        """

        self.repetition_period = repetition_period
        self.init_state = init_state
        self.external_error_probs = external_error_probs

        # Apply time window selection if specified
        if center_time is not None and window is not None:
            # Create time array for original measurements
            times = np.array([i * repetition_period for i in range(len(measurements))])

            t_start, t_end = center_time + window[0], center_time + window[1]
            
            # Create selection mask
            sel = (times >= t_start) & (times <= t_end)
            
            # Apply selection to measurements
            self.measurements = [measurements[i] for i in range(len(measurements)) if sel[i]]
            self.time_offset = np.where(sel)[0][0] if np.any(sel) else 0  # Store offset for time calculation
        else:
            self.measurements = measurements
            self.time_offset = 0

        self.errors = {} 
        self.analysis_type = None  # Store the type of analysis used

        self.threshold = threshold  # threshold to determine high or low period

        self.verbose = verbose  # verbosity level
    
    def analyze_errors(self):
        '''
        Method to compute the number of errors according to each case.
        '''
        if self.repetition_period >= self.threshold:
            self.analysis_type = "high_period"
            self._analyze_high_period()
        else:
            self.analysis_type = "low_period"
            self._analyze_low_period()




    def _analyze_high_period(self):

        '''
        Method to compute the number of errors when there is enough time
        between measurements for the state to go back to |0>, i.e.,
        when repetition_period >= 150 mus. 

        Theoretical histograms are only computed if the number of qubits is 
        less than 7.
        '''
            

        expected_value = '0' if self.init_state == 0 else '1'
        num_measurements = len(self.measurements)
        num_qubits = len(self.measurements[0])

        # Join all bitstrings and convert ASCII '0'/'1' -> integers 0/1
        flat = np.frombuffer(''.join(self.measurements).encode('ascii'), dtype=np.uint8) - ord('0')

        # Reshape into (n_measurements, n_qubits)
        meas_array = flat.reshape(num_measurements, num_qubits)

        # Expected array: all 0s or all 1s
        expected_bit = 0 if expected_value == '0' else 1
        expected_array = np.full_like(meas_array, expected_bit, dtype=np.uint8)

        # Compute mismatch mask (1 if error, 0 if correct)
        error_mask = (meas_array != expected_array).astype(np.uint8)

        # Compute total errors per measurement (sum over qubits)
        total_errors = error_mask.sum(axis=1)

        # Compute errors per qubit (sum over measurements)
        errors_per_qubit = error_mask.sum(axis=0)

        # Compute times
        times = (np.arange(num_measurements) + self.time_offset) * self.repetition_period

        # Generate histogram of error patterns (only if small enough)
        error_pattern_histogram = Counter()
        if num_qubits <= 6:
            for mask in error_mask:
                errored_qubits = tuple(np.nonzero(mask)[0])
                error_pattern_histogram[errored_qubits] += 1
        else:
            error_pattern_histogram = None  # Skip to save time

        # Build output
        self.errors["high_period_errors"] = {
            "total_errors": total_errors.tolist(),
            "times": times.tolist(),
            "errors_per_qubit": errors_per_qubit.tolist(),
            "states": self.measurements
        }

        # If feasible, compute theoretical histogram
        if num_qubits <= 6:
            total_counts = num_measurements
            p_qubit_error = self._get_error_probabilities(errors_per_qubit, total_counts)

            theoretical_histogram = {}
            for r in range(num_qubits + 1):
                for qubit_combo in combinations(range(num_qubits), r):
                    prob = np.prod([
                        p_qubit_error[q] if q in qubit_combo else (1 - p_qubit_error[q])
                        for q in range(num_qubits)
                    ])
                    theoretical_histogram[qubit_combo] = prob * total_counts

            self.errors["high_period_errors"]["error_patterns"] = dict(error_pattern_histogram)
            self.errors["high_period_errors"]["theoretical_histogram"] = theoretical_histogram

        else:
            # If skipped due to high qubit count
            self.errors["high_period_errors"]["error_patterns"] = None
            self.errors["high_period_errors"]["theoretical_histogram"] = None

        self.total_counts = num_measurements
    
    def _analyze_low_period(self):
        '''
        Method to compute the number of errors when X gates are applied 
        continuously, i.e., when repetition_period < 500 mus.

        The number of errors is defined from the previous measurement.

        A classification of the errors is made here in True and False errors,
        according to the already given definitions.

        Theoretical histograms are only computed if the number of qubits is
        less than 7.
        '''

        num_measurements = len(self.measurements)
        num_qubits = len(self.measurements[0])

        # Convert measurements to NumPy array
        meas_array = np.frombuffer(''.join(self.measurements).encode('ascii'), dtype=np.uint8) - ord('0')
        meas_array = meas_array.reshape(num_measurements, num_qubits)

        # Compute expected measurements
        expected_strings = self._compute_expected_sequences()
        expected_array = np.frombuffer(''.join(expected_strings).encode('ascii'), dtype=np.uint8) - ord('0')
        expected_array = expected_array.reshape(num_measurements - 1, num_qubits)

        # Compute masks for true/false errors
        true_errors_mask = (expected_array == 0) & (meas_array[1:] == 1)
        false_errors_mask = (expected_array == 1) & (meas_array[1:] == 0)

        # Total errors per measurement
        true_errors = true_errors_mask.sum(axis=1)
        false_errors = false_errors_mask.sum(axis=1)
        total_errors = true_errors + false_errors

        # Errors per qubit
        errors_per_qubit = (true_errors_mask | false_errors_mask).sum(axis=0)

        # Compute times
        times = (np.arange(num_measurements - 1) + self.time_offset) * self.repetition_period

        # Generate histogram of error patterns (only if small number of qubits)
        error_pattern_histogram = None
        if num_qubits <= 6:
            error_pattern_histogram = Counter()
            for mask in (true_errors_mask | false_errors_mask):
                error_pattern_histogram[tuple(np.nonzero(mask)[0])] += 1

        # Compute theoretical histogram only if num_qubits <= 6
        theoretical_histogram = None
        if num_qubits <= 6:
            total_counts = num_measurements - 1
            p_qubit_error = self._get_error_probabilities(errors_per_qubit, total_counts)
            theoretical_histogram = {}
            for r in range(num_qubits + 1):
                for qubit_combo in combinations(range(num_qubits), r):
                    prob = np.prod([
                        p_qubit_error[q] if q in qubit_combo else (1 - p_qubit_error[q])
                        for q in range(num_qubits)
                    ])
                    theoretical_histogram[qubit_combo] = prob * total_counts

        # Store results
        self.errors["low_period_errors"] = {
            "true_errors": true_errors.tolist(),
            "false_errors": false_errors.tolist(),
            "total_errors": total_errors.tolist(),
            "times": times.tolist(),
            "errors_per_qubit": errors_per_qubit.tolist(),
            "states": self.measurements[1:],
            "expected": expected_strings,
            "error_patterns": dict(error_pattern_histogram) if error_pattern_histogram else None,
            "theoretical_histogram": theoretical_histogram
        }

        self.total_counts = num_measurements - 1


    
    def _get_error_probabilities(self, errors_per_qubit, total_counts):
        """
        Get error probabilities either from external source or computed from data.
        
        Parameters:
        - errors_per_qubit: array of error counts per qubit
        - total_counts: total number of measurements
        
        Returns:
        - p_qubit_error: array of error probabilities for each qubit
        """
        if self.external_error_probs is not None:
            # Use external probabilities
            p_qubit_error = np.array(self.external_error_probs)
            
            # Validate dimensions
            if len(p_qubit_error) != len(errors_per_qubit):
                raise ValueError(f"External error probabilities length ({len(p_qubit_error)}) "
                            f"must match number of qubits ({len(errors_per_qubit)})")
            
            # Validate probability range
            if np.any(p_qubit_error < 0) or np.any(p_qubit_error > 1):
                raise ValueError("External error probabilities must be between 0 and 1")
                
            return p_qubit_error
        else:
            # Compute from observed data
            return np.array(errors_per_qubit) / total_counts
    
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


    def get_qubit_error_probabilities(self, center_time, window=(-0.01, 0.015), bin_size=4, qubit: int = None):
        """
        Compute error probabilities per qubit in a given time window, binned over a fixed number of measurements.

        This version supports both high_period and low_period analysis modes:

        - high_period: compares each measurement to the fixed init_state.
        - low_period: compares each measurement (from index 1) to the expected sequence
                      computed from the previous measurement (same logic as _analyze_low_period).

        Parameters:
        - center_time: float, center of the time window (absolute time, seconds)
        - window: tuple (t_before, t_after) relative to center_time
        - bin_size: int, number of measurements per bin (if selected measurements < bin_size we return a single bin)
        - qubit: optional int, index of a single qubit to compute. If None, compute for all qubits.

        Returns:
            times_bin: np.ndarray of time per bin (center)
            p_qubits_bin: np.ndarray (n_bins, n_qubits) of average error probability per qubit
                          or np.ndarray (n_bins,) if qubit is specified
        """
        if self.analysis_type is None:
            raise RuntimeError("Run analyze_errors() first to determine analysis_type.")

        # Convert measurements to array
        num_measurements = len(self.measurements)
        if num_measurements == 0:
            return np.array([]), np.array([])

        n_qubits = len(self.measurements[0])
        meas_array = np.frombuffer(''.join(self.measurements).encode('ascii'), dtype=np.uint8) - ord('0')
        meas_array = meas_array.reshape(num_measurements, n_qubits)

        if self.analysis_type == "high_period":
            # expected bit for high_period mode must be provided
            if self.init_state is None:
                raise RuntimeError("init_state must be set for high_period analysis to compute error probabilities.")
            expected_bit = 0 if self.init_state == 0 else 1
            expected_array = np.full_like(meas_array, expected_bit, dtype=np.uint8)

            # Compute error mask (1 if error, 0 if correct)
            error_mask = (meas_array != expected_array).astype(np.uint8)

            # Absolute times for the measurements (accounting for any time_offset from constructor selection)
            times = (np.arange(num_measurements) + self.time_offset) * self.repetition_period

        elif self.analysis_type == "low_period":
            # low_period: expected sequence depends on previous measurement -> compares meas_array[1:] with expected_array
            if num_measurements < 2:
                return np.array([]), np.array([])

            expected_strings = self._compute_expected_sequences()  # length num_measurements-1
            expected_array = np.frombuffer(''.join(expected_strings).encode('ascii'), dtype=np.uint8) - ord('0')
            expected_array = expected_array.reshape(num_measurements - 1, n_qubits)

            # observed for comparison are meas_array[1:]
            observed = meas_array[1:, :]

            # true/false error masks as in _analyze_low_period
            true_errors_mask = (expected_array == 0) & (observed == 1)
            false_errors_mask = (expected_array == 1) & (observed == 0)

            error_mask = (true_errors_mask | false_errors_mask).astype(np.uint8)

            # times correspond to measurements 1..end (shifted)
            times = (np.arange(num_measurements - 1) + self.time_offset) * self.repetition_period

        else:
            raise RuntimeError(f"Unknown analysis_type: {self.analysis_type}")

        # Select measurements in requested window
        tmin = center_time + window[0]
        tmax = center_time + window[1]
        mask = (times >= tmin) & (times <= tmax)

        selected_count = np.count_nonzero(mask)
        print("Selected measurements:", selected_count)

        if selected_count == 0:
            return np.array([]), np.array([])

        selected_errors = error_mask[mask]
        selected_times = times[mask]

        # If user requested a single qubit, reduce shape
        if qubit is not None:
            if qubit < 0 or qubit >= n_qubits:
                raise IndexError("qubit index out of range")
            selected_errors_single = selected_errors[:, qubit]
            if selected_errors_single.shape[0] < bin_size:
                n_bins = 1
            else:
                n_bins = selected_errors_single.shape[0] // bin_size
            times_bin = np.zeros(n_bins)
            p_bin = np.zeros(n_bins)
            for i in range(n_bins):
                start = i * bin_size
                end = start + bin_size if n_bins > 1 else selected_errors_single.shape[0]
                block = selected_errors_single[start:end]
                p_bin[i] = block.mean()
                times_bin[i] = selected_times[start:end].mean()
            return times_bin, p_bin

        # Multi-qubit case
        if selected_errors.shape[0] < bin_size:
            n_bins = 1
        else:
            n_bins = selected_errors.shape[0] // bin_size

        p_qubits_bin = np.zeros((n_bins, n_qubits))
        times_bin = np.zeros(n_bins)
        for i in range(n_bins):
            start = i * bin_size
            end = start + bin_size if n_bins > 1 else selected_errors.shape[0]
            block = selected_errors[start:end, :]
            p_qubits_bin[i, :] = block.mean(axis=0)
            times_bin[i] = selected_times[start:end].mean()

        return times_bin, p_qubits_bin

    
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
            print(f"(repetition_period ≥ {self.threshold*1e-6} μs). The error dictionary contains:\n")
            print(f"- Key: 'high_period_errors' → Returns two arrays with {len(self.errors['high_period_errors']['total_errors'])} measurements:")
            print("  - 'error_counts': Array with the number of errors associated with each measurement.")
            print("  - 'times': Array with the times of each measurement (calculated as i * repetition_period).")
            print("\n")
        else:
            print("The analysis was performed using the *low repetition period* approach ")
            print(f"(repetition_period < {self.threshold*1e-6} μs). The error dictionary contains:\n")
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
        



