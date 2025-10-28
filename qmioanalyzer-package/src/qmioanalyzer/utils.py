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

            # Skip header lines - handle both space and tab separated
            line_lower = line.strip().lower()
            if (line_lower.startswith(("batch", "repetition", "state")) or 
                "state" in line_lower and "abstime" in line_lower):
                continue  

            # Parse data lines - handle both space and tab separation
            parts = line.strip().split()  # This handles both spaces and tabs
            
            if len(parts) == 2:
                # Format: State AbsTime[s] (tab or space separated)
                try:
                    bits.append(parts[0])  # State is first column
                    abstime.append(float(parts[1]))  # AbsTime is second column
                except Exception:
                    # fallback: treat as old format [index, bitstring]
                    bits.append(parts[1])
            elif len(parts) >= 3:
                # Format: Batch Shot State AbsTime[s] or other multi-column
                try:
                    batch.append(int(parts[0]))
                    shot.append(int(parts[1]))
                    bits.append(parts[2])
                    if len(parts) > 3:
                        abstime.append(float(parts[3]))
                except Exception:
                    # fallback: assume state is in appropriate column
                    if len(parts) == 3:
                        # Could be: State AbsTime[s] with extra column
                        bits.append(parts[0])
                        try:
                            abstime.append(float(parts[1]))
                        except:
                            pass
                    else:
                        bits.append(parts[2])  # default to third column
            elif len(parts) == 1:
                # Single column - just the state
                bits.append(parts[0])

    extra = {}
    if batch:   extra["batch"] = np.array(batch)
    if shot:    extra["shot"] = np.array(shot)
    if abstime: extra["abstime"] = np.array(abstime)

    return np.array(bits), repetition_period, init_state, backend, extra