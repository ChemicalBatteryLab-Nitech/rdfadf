import numpy as np
from ase.neighborlist import NeighborList
from scipy.ndimage import gaussian_filter1d
from .label_utils import interpret_label


def compute_general_rdf(
    atoms,
    center_label,
    neighbor_label,
    rcut=10.0,
    bins=200,
    sigma=2.0
):
    """
    RDF with:
        - correct PBC handling
        - center-atom normalization
        - cumulative counts (NOT integral)
    """

    atoms.pbc = True

    positions = atoms.positions
    cell = atoms.cell
    symbols = atoms.get_chemical_symbols()
    n = len(atoms)

    # label expansion
    center_set = interpret_label(center_label, symbols)
    neighbor_set = interpret_label(neighbor_label, symbols)

    # center atoms
    center_indices = [i for i, s in enumerate(symbols) if s in center_set]
    n_center = len(center_indices)
    if n_center == 0:
        raise ValueError("No atoms match center_label")

    # neighbor list (PBC)
    cutoffs = [rcut] * n
    nl = NeighborList(
        cutoffs,
        self_interaction=False,
        bothways=True,
        skin=0.0
    )
    nl.update(atoms)

    # distance collection
    dist_list = []
    for i in center_indices:
        neigh, offsets = nl.get_neighbors(i)
        for j, off in zip(neigh, offsets):
            if symbols[j] not in neighbor_set:
                continue
            dr = positions[j] + off @ cell - positions[i]
            r = np.linalg.norm(dr)
            if r < rcut:
                dist_list.append(r)

    # histogram
    hist, edges = np.histogram(dist_list, bins=bins, range=(0, rcut))
    r = 0.5 * (edges[:-1] + edges[1:])
    g_raw = hist.astype(float)

    # normalize by number of center atoms
    g_raw_norm = g_raw / n_center

    # smoothing
    g_smooth = gaussian_filter1d(g_raw_norm, sigma=sigma)

    # 4*pi*r^2 normalization
    eps = 1e-8
    denom = 4.0 * np.pi * (r + eps) ** 2
    g_raw2 = g_raw_norm / denom
    g_smooth2 = g_smooth / denom

    # ===== cumulative counts (NOT integral!) =====
    cumulative_raw = np.cumsum(g_raw_norm)
    cumulative_smooth = np.cumsum(g_smooth)
    #cumulative_raw2 = np.cumsum(g_raw2)
    #cumulative_smooth2 = np.cumsum(g_smooth2)

    return (
        r,
        g_raw_norm,
        g_smooth,
        g_raw2,
        g_smooth2,
        cumulative_raw,
        cumulative_smooth,
        #cumulative_raw2,
        #cumulative_smooth2
    )