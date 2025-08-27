import numpy as np


def read_data(filename):
    '''
    Reads a data file with optional header lines beginning with '#'.

    Returns:
    - bits (np.array of str): list of bitstrings ("State" column)
    - repetition_period (float or None)
    - init_state (int, None, or str if unrecognized)
    - backend (str or None): "qmio", "fakeqmio", or None if unknown
    - extra (dict): may contain keys "batch", "shot", "abstime"
                    each as np.array, or empty if not available
    '''
    bits = []
    repetition_period = None
    init_state = None
    backend = None

    # new arrays for optional columns
    batch = []
    shot = []
    abstime = []

    with open(filename, "r") as file:
        for line in file:
            if line.startswith("#"):  # metadata
                line = line.strip()
                if "repetition_period" in line:
                    try:
                        repetition_period = float(line.split(":")[1].strip().split()[0])
                        backend = "qmio"
                    except Exception:
                        pass
                elif "delay" in line and "circuit_time" in line:  # fakeqmio case
                    try:
                        parts = line.split(",")
                        circuit_time_str = parts[1].split(":")[1].strip().split()[0]
                        repetition_period = float(circuit_time_str)
                        backend = "fakeqmio"
                    except Exception:
                        pass
                elif "state initialized to" in line:
                    val = line.split(":")[1].strip()
                    if val.lower() == "none":
                        init_state = None
                    else:
                        try:
                            init_state = int(val)
                        except Exception:
                            init_state = val
                continue

            # skip header line ("Batch Shot State AbsTime[s]")
            if line.strip().lower().startswith(("batch", "repetition", "state")):
                continue  

            parts = line.strip().split()
            if len(parts) == 2:
                # old format: [index, bitstring]
                bits.append(parts[1])
            elif len(parts) >= 3:
                # new format: Batch Shot State AbsTime[s]
                try:
                    batch.append(int(parts[0]))
                    shot.append(int(parts[1]))
                    bits.append(parts[2])
                    if len(parts) > 3:
                        abstime.append(float(parts[3]))
                except Exception:
                    # fallback: only bits
                    bits.append(parts[2])

    extra = {}
    if batch:   extra["batch"] = np.array(batch)
    if shot:    extra["shot"] = np.array(shot)
    if abstime: extra["abstime"] = np.array(abstime)

    return np.array(bits), repetition_period, init_state, backend, extra


