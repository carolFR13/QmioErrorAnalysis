from .error_analyzer import ErrorAnalyzer
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from typing import Literal

class ErrorVisualizer:
    ''''''

    def __init__(self, error_obj, error_obj_2 = None, normalize = False):

        self.error_obj = error_obj
        self.init_state = error_obj.init_state
        self.repetition_period = error_obj.repetition_period
        self.error_dict = error_obj.get_errors()


        self.subplot = True if error_obj_2 else False
        self.error_obj_2 = error_obj_2 if error_obj_2 else None
        self.init_state_2 = error_obj_2.init_state if error_obj_2 else None
        self.repetition_period_2 = error_obj_2.repetition_period if error_obj_2 else None
        self.error_dict_2 = error_obj_2.get_errors() if error_obj_2 else None


        self.normalize = normalize

    def _get_error_data(self, period_type: str, obj: Literal[0, 1] = 0) -> dict:
        """
        Returns the relevant error dictionary depending on repetition period and object.

        period_type: "low" or "high"
        obj: 0 or 1 -> select which ErrorAnalyzer object to return data from
        """
        error_dict = self.error_dict if obj == 0 else self.error_dict_2

        if period_type == "high":
            return error_dict['high_period_errors']
        elif period_type == "low":
            return error_dict['low_period_errors']
        else:
            raise ValueError("Unknown period_type. Use 'low' or 'high'.")

    def get_states_histogram(self, period: str, obj: int = 0):
        '''
        method to obtain the preferred states of the qubits
        '''
        data = self._get_error_data(period, obj)
        states = data['states']
        total_errors = data['total_errors']

        # defining a dictionary to count the number of errors assigned to each state
        state_error_counts = defaultdict(int) 

        # assigning each state with the corresponding number of errors
        for state, error in zip(states, total_errors):
            state_error_counts[state] += error

        # we order the states according to their total number of errors
        states_sorted = sorted(state_error_counts.keys(), key=lambda x: state_error_counts[x], reverse=True)
        errors_sorted = [state_error_counts[state] for state in states_sorted]

        return states_sorted, errors_sorted

    def get_error_patterns(self, period: str = None, obj: int = 0, type: str = "experimental", threshold: float = 0.0):
        '''
        type can be "theoretical", "experimental" or "both".
        '''

        data = self._get_error_data(period, obj)

        if type == "theoretical":
            patterns = data['theoretical_histogram']
        elif type == "experimental":
            patterns = data['error_patterns']
        else:
            patterns_exp = data['error_patterns']
            patterns_theo = data['theoretical_histogram']

            all_keys = set(patterns_exp.keys()).union(patterns_theo.keys())
            total_exp = sum(patterns_exp.values())
            total_theo = sum(patterns_theo.values())

            print(f"[DEBUG] Experimental keys: {len(patterns_exp.keys())}")
            print(f"[DEBUG] Theoretical keys: {len(patterns_theo.keys())}")
            print(f"[DEBUG] All keys: {len(all_keys)}")
            print(f"[DEBUG] Total exp: {total_exp}, Total theo: {total_theo}")
            
            # Filter based on EITHER experimental OR theoretical significance
            filtered_keys = []
            for k in all_keys:
                exp_freq = patterns_exp.get(k, 0) / total_exp if total_exp > 0 else 0
                theo_freq = patterns_theo.get(k, 0) / total_theo if total_theo > 0 else 0
                
                # Include if either experimental OR theoretical exceeds threshold
                if exp_freq >= threshold or theo_freq >= threshold:
                    filtered_keys.append(k)

            filtered_keys = sorted(filtered_keys, key=lambda k: patterns_exp.get(k, 0), reverse=True)

            print(f"[INFO] Bins before cut: {len(all_keys)}")
            print(f"[INFO] Bins after cut: {len(filtered_keys)}")

            labels_both = [str(list(k)) for k in filtered_keys]
            observed_counts = [patterns_exp.get(k, 0) for k in filtered_keys]
            theoretical_counts = [patterns_theo.get(k, 0) for k in filtered_keys]
            
            return labels_both, observed_counts, theoretical_counts

        total = sum(patterns.values())
        filtered = [(k, v) for k, v in patterns.items() if v / total >= threshold]
        filtered_sorted = sorted(filtered, key=lambda x: x[1], reverse=True)

        print(f"[INFO] Bins before cut: {len(patterns)}")
        print(f"[INFO] Bins after cut: {len(filtered_sorted)}")

        labels = [str(list(k)) for k, _ in filtered_sorted]
        freqs = [v for _, v in filtered_sorted]
        return labels, freqs

    def get_error_counts(self, period: str, obj: int = 0):

        if period == "high":
            data = self._get_error_data("high", obj)
            total_errors = data['total_errors']
            return total_errors

        elif period == "low":
            data = self._get_error_data("low", obj)
            total_errors = data['total_errors']
            true_errors = data['true_errors']
            false_errors = data['false_errors']

            return total_errors, true_errors, false_errors

    def get_error_counts_per_qubit(self, period: str, obj: int = 0):

        data = self._get_error_data(period, obj)
        errors_per_qubit = data['errors_per_qubit']

        return errors_per_qubit

    @staticmethod
    def plot_histogram(labels, values, ax=None, title="", color="blue", rotation=90):
        '''
        Generic method to plot a histogram of errors.
        Appropiate for states plot and error patterns plot.
        '''
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.bar(labels, values, color=color)
        ax.set_title(title)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=rotation)
        return ax



