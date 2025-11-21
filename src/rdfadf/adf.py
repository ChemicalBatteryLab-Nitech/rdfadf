import numpy as np
from ase.neighborlist import NeighborList
from scipy.ndimage import gaussian_filter1d
from .label_utils import interpret_label

def compute_general_adf(
    atoms,
    center_label="All",
    neigh1_label="All",
    neigh2_label="All",
    rcut=6.0,
    bins=180,
    sigma=None,
    r_min=0.0,
    r_max=None,
):
    """
    Angle Distribution Function (ADF) with:
        - PBC support
        - element/block/group/period filtering
        - normalization by number of center atoms
        - radial window [r_min, r_max] for neighbors j,k

    Parameters
    ----------
    atoms : ASE Atoms
    center_label : str
        Classification for central atom (element, All, s/p/d/f, G#, P#)
    neigh1_label : str
        Classification for neighbor j
    neigh2_label : str
        Classification for neighbor k
    rcut : float
        Cutoff for neighbor list construction
    bins : int
        Number of angle bins between 0 and pi
    sigma : float or None
        Gaussian smoothing width (None -> no smoothing)
    r_min : float
        Minimum distance from center to neighbors j,k to be included
    r_max : float or None
        Maximum distance from center to neighbors j,k to be included
        If None, r_max = rcut

    Returns
    -------
    theta : (bins,) array
        Angle bin centers (radians)
    adf_raw_norm : (bins,) array
        ADF histogram per center atom
    adf_smooth : (bins,) array
        Smoothed ADF histogram per center atom
    """

    atoms.pbc = True

    positions = atoms.positions
    cell = atoms.cell
    symbols = atoms.get_chemical_symbols()
    n = len(atoms)

    if r_max is None:
        r_max = rcut

    # expand labels
    center_set = interpret_label(center_label, symbols)
    neigh1_set = interpret_label(neigh1_label, symbols)
    neigh2_set = interpret_label(neigh2_label, symbols)

    # center atoms
    center_indices = [i for i, s in enumerate(symbols) if s in center_set]
    n_center = len(center_indices)
    if n_center == 0:
        raise ValueError("No atoms match center_label")

    # neighbor list (PBC-aware)
    cutoffs = [rcut] * n
    nl = NeighborList(
        cutoffs,
        self_interaction=False,
        bothways=True,
        skin=0.0,
    )
    nl.update(atoms)

    angles = []

    # angle loop
    for i in center_indices:
        neigh, offs = nl.get_neighbors(i)
        m = len(neigh)

        for a in range(m):
            j = neigh[a]
            if symbols[j] not in neigh1_set:
                continue

            dr1 = positions[j] + offs[a] @ cell - positions[i]
            d1 = np.linalg.norm(dr1)
            if not (r_min <= d1 <= r_max):
                continue

            for b in range(a + 1, m):
                k = neigh[b]
                if symbols[k] not in neigh2_set:
                    continue

                dr2 = positions[k] + offs[b] @ cell - positions[i]
                d2 = np.linalg.norm(dr2)
                if not (r_min <= d2 <= r_max):
                    continue

                # angle at center: j - i - k
                cos_th = np.dot(dr1, dr2) / (
                    np.linalg.norm(dr1) * np.linalg.norm(dr2)
                )
                cos_th = np.clip(cos_th, -1.0, 1.0)
                theta = np.arccos(cos_th)

                angles.append(theta)

    # histogram
    hist, edges = np.histogram(angles, bins=bins, range=(0, np.pi))
    theta = 0.5 * (edges[:-1] + edges[1:])
    adf_raw = hist.astype(float)

    # normalize per center atom
    adf_raw_norm = adf_raw / n_center

    # smoothing
    if sigma is not None:
        adf_smooth = gaussian_filter1d(adf_raw_norm, sigma=sigma)
    else:
        adf_smooth = adf_raw_norm.copy()

    return theta, adf_raw_norm, adf_smooth
