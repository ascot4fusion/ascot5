import numpy as np

def loads(ids, energy, weight, area):
    """
    Get wall loads per triangle.

    Returns:
            ids
            edepo
            eload
            pdepo
            pload
            mdepo
            iangle
    """
    # Indices of wetted triangles
    wetted = np.unique(ids)
    area   = area[wetted-1]
    edepo  = np.zeros(wetted.shape)
    eload  = np.zeros(wetted.shape)
    pdepo  = np.zeros(wetted.shape)
    eload  = np.zeros(wetted.shape)
    mdepo  = np.zeros(wetted.shape)
    iangle = np.zeros(wetted.shape)

    # Sort markers by wall tile ID to make processing faster
    idx    = np.argsort(ids)
    ids    = ids[idx]
    energy = energy[idx]
    weight = weight[idx]

    idx = np.append(np.argwhere(ids[1:] - ids[:-1]).ravel(), ids.size-1)

    i0 = 0
    for i in range(wetted.size):
        i1 = idx[i] + 1
        mdepo[i]  = i1 - i0
        pdepo[i]  = np.sum(weight[i0:i1])
        edepo[i]  = np.sum(energy[i0:i1]*weight[i0:i1])
        iangle[i] = 0 #TBD
        i0 = i1
    pload = pdepo / area
    eload = edepo / area

    return (wetted, edepo, eload, pdepo, pload, mdepo, iangle)


def figuresofmerit(ids, energy, weight, area, plim=None):
    """
    Return figures of merit (0D quantities) for wall loads.

    Returns:
    peak power/energy load
    peak particle load
    number of markers
    total wetted area
    wetted area above limi
    """
    wetted, _, eload, _, pload, mdepo, _ = loads(ids, energy, weight, area)

    wetted_total = np.sum(area[wetted-1])
    energy_peak  = np.amax(eload)

    return (wetted_total, energy_peak)


def loadvsarea(ids, energy, weight, area):
    """
    Return a curve showing affected area as a function of (minimum) load.
    """
    wetted, _, eload, _, _, _, _ = loads(ids, energy, weight, area)
    idx = np.argsort(-eload)
    y = (eload[idx])
    x = np.cumsum(w[idx])
