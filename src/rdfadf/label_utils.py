from .periodic_info import periodic_info


def element_to_label(symbol, mode):
    block, group, period = periodic_info[symbol]

    if mode == "element":
        return symbol
    elif mode == "block":
        return block
    elif mode == "group":
        return f"G{group}"
    elif mode == "period":
        return f"P{period}"
    else:
        raise ValueError("Invalid mode")

def interpret_label(label, symbols):
    """
    Expand an input label into a set of element symbols contained in `atoms`.

    Supported labels
    ----------------
    - "All"              : all element symbols present in `symbols`
    - Element symbol     : e.g., "Li", "O", "Co"
    - Block label        : "s", "p", "d", "f"
    - Group label        : "G#" (e.g., "G1", "G16")
    - Period label       : "P#" (e.g., "P2", "P4")

    Parameters
    ----------
    label : str
        Classification label for selecting elements.
    symbols : list of str
        List of element symbols from the Atoms object.

    Returns
    -------
    set of str
        Set of element symbols that match the given label.
    """
    unique_syms = set(symbols)
    valid = set(periodic_info.keys())

    # 1) All elements present in the structure
    if label == "All":
        return unique_syms

    # 2) Exact element symbol
    if label in valid:
        return {label}

    # 3) Block: s, p, d, f
    if label in ["s", "p", "d", "f"]:
        return {s for s in unique_syms if periodic_info[s][0] == label}

    # 4) Group: G1, G2, ...
    if label.startswith("G"):
        g = int(label[1:])
        return {s for s in unique_syms if periodic_info[s][1] == g}

    # 5) Period: P1, P2, ...
    if label.startswith("P"):
        p = int(label[1:])
        return {s for s in unique_syms if periodic_info[s][2] == p}

    raise ValueError(f"Unknown label: {label}")