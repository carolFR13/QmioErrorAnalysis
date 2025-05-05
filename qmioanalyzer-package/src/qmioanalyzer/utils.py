import numpy as np

def read_data(filename):
    '''
    Function to read a file of data with two separated columns.

    It returns a numpy array with the values of the second column as 
    strings.
    '''
    
    bits = []
    with open(filename, "r") as file:
        next(file) # jumping to the second line. first is header
        for line in file:
            parts = line.split()
            if len(parts) == 2:
                bits.append(parts[1])  # we save the bits string in type string

    return np.array(bits)  