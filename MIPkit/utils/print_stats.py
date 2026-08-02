import numpy as np
import os
from scipy.spatial.distance import cdist


def print_update(i, nmono=0, npoly=0, mws=0, init=False):
    if init:
        a = "Step".ljust(6)
        b = "Nmono".rjust(6)
        c = "Npoly".rjust(6)
        d1 = "MW-Navg".rjust(12)
        d2 = "MW-Wavg".rjust(12)
        d3 = "MWmax".rjust(12)
        d4 = "MWmin".rjust(12)
        print("======================================================================")
        print("  %s%s%s%s%s%s%s" % (a, b, c, d1, d2, d3, d4))
        print("======================================================================")
    else:
        a = str(i).ljust(6)
        b = str(nmono).rjust(6)
        c = str(npoly).rjust(6)
        d1 = str(round(np.mean(mws), 2)).rjust(12)

        mww = 0
        for mw in mws:
            mww += mw / np.sum(mws) * mw

        d2 = str(round(mww, 2)).rjust(12)
        d3 = str(round(np.max(mws), 2)).rjust(12)
        d4 = str(round(np.min(mws), 2)).rjust(12)
        print("  %s%s%s%s%s%s%s" % (a, b, c, d1, d2, d3, d4))

    return


def print_stats(workdir, i=0, nmono=0, npoly=0, mws=0, init=False, precomplex="", template=""):
    if init:
        stats = open(f"{workdir}/stats_polymerization.txt", "w+")
        a = "Step".ljust(6)
        b = "Nmono".rjust(6)
        c = "Npoly".rjust(6)
        d1 = "MW-Navg".rjust(12)
        d2 = "MW-Wavg".rjust(12)
        d3 = "MWmax".rjust(12)
        d4 = "MWmin".rjust(12)

        stats.write(f"Precomplex : {precomplex}\nTemplate : {template}\n")

    else:
        stats = open(f"{workdir}/stats_polymerization.txt", "a")
        a = str(i).ljust(6)
        b = str(nmono).rjust(6)
        c = str(npoly).rjust(6)
        d1 = str(round(np.mean(mws), 2)).rjust(12)

        mww = 0
        for mw in mws:
            mww += mw / np.sum(mws) * mw

        d2 = str(round(mww, 2)).rjust(12)
        d3 = str(round(np.max(mws), 2)).rjust(12)
        d4 = str(round(np.min(mws), 2)).rjust(12)

    stats.write("%s%s%s%s%s%s%s\n" % (a, b, c, d1, d2, d3, d4))
    stats.close()

    return


def print_mw_distributions(
    workdir, i=0, totmon=0, mws=0, init=False, precomplex="", template=""
):
    if init:
        stats = open(f"{workdir}/stats_mw_distributions.txt", "w+")
        a = "Step".ljust(6)
        b = "MW...".rjust(12)

        stats.write(f"Precomplex : {precomplex}\nTemplate : {template}\n")

    else:
        stats = open(f"{workdir}/stats_mw_distributions.txt", "a")

        allmw = np.zeros(totmon)

        for j, mw in enumerate(mws):
            allmw[j] = mw

        a = str(i).ljust(6)

        b = ""
        for mw in allmw:
            b += str(round(mw, 2)).rjust(12)

    stats.write("%s%s\n" % (a, b))
    stats.close()

    return


def print_bondstats(
    workdir,
    statname,
    step=0,
    species_count=0,
    fma=0,
    positionsa=0,
    atomsa=0,
    fmb=0,
    positionsb=0,
    atomsb=0,
    prot_coords=0,
    prot_res=0,
    prot_atom=0,
    init=False,
    precomplex="",
    template="",
):
    if init:
        poly_stats = open(os.path.join(workdir, statname), "w+")
        a = "Step".ljust(6)
        b = "RXN".ljust(6)
        c0 = "Dist".ljust(9)
        c1 = "FM".ljust(9)
        c2 = "Atom".ljust(9)
        c3 = "Res".ljust(9)
        c4 = "Atom".ljust(9)

        poly_stats.write(f"Precomplex : {precomplex}\nTemplate : {template}\n")

        poly_stats.write("%s%s%s%s%s%s%s\n" % (a, b, c0, c1, c2, c3, c4))
        poly_stats.close()

    else:
        distsa = cdist(positionsa, prot_coords)
        distsb = cdist(positionsb, prot_coords)

        loc = []

        out_dist = 0
        cfm = "NaN"
        catoms = []

        if distsa.min() < distsb.min():
            out_dist = distsa.min()
            for j in np.where(distsa == distsa.min()):
                loc.append(int(j[0]))
            cfm = fma
            catoms = atomsa

        else:
            out_dist = distsb.min()
            for j in np.where(distsb == distsb.min()):
                loc.append(int(j[0]))
            cfm = fmb
            catoms = atomsb

        poly_stats = open(os.path.join(workdir, statname), "a+")
        a = str(step).ljust(6)
        b = str(species_count).ljust(6)
        c0 = str(round(out_dist, 2)).ljust(9)
        c1 = cfm.ljust(9)
        c2 = (catoms[loc[0]] + str(loc[0])).ljust(9)
        c3 = prot_res[loc[1]].ljust(9)
        c4 = (prot_atom[loc[1]] + str(loc[1])).ljust(9)
        poly_stats.write("%s%s%s%s%s%s%s\n" % (a, b, c0, c1, c2, c3, c4))
        poly_stats.close()

    return
