import argparse
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial.distance import cdist
import pandas as pd

from MIPkit.constants.constants import encode_fms
from MIPkit.utils.wildcard import wildcard
from MIPkit.utils.utils import generate_recipe
from MIPkit.utils.read_yaml import read_yaml
from MIPkit.utils.read_pdb import read_pdb
from MIPkit.utils.parser import sanitize_visual


def add_headers(
        fig,
        *,
        row_headers=None,
        col_headers=None,
        row_pad=1,
        col_pad=5,
        rotate_row_headers=True,
        **text_kwargs
):
    # Based on https://stackoverflow.com/a/25814386

    axes = fig.get_axes()

    for ax in axes:
        sbs = ax.get_subplotspec()

        # Putting headers on cols
        if (col_headers is not None) and sbs.is_first_row():
            ax.annotate(
                col_headers[sbs.colspan.start],
                xy=(0.5, 1),
                xytext=(0, col_pad),
                xycoords="axes fraction",
                textcoords="offset points",
                ha="center",
                va="baseline",
                **text_kwargs,
            )

        # Putting headers on rows
        if (row_headers is not None) and sbs.is_first_col():
            ax.annotate(
                row_headers[sbs.rowspan.start],
                xy=(0, 0.5),
                xytext=(-ax.yaxis.labelpad - row_pad, 0),
                xycoords=ax.yaxis.label,
                textcoords="offset points",
                ha="right",
                va="center",
                rotation=rotate_row_headers * 90,
                **text_kwargs,
            )

def trajectory_distances(templatename=None,cplxname=None,recipe=None):

    enfm = encode_fms()

    # if args is not None and args.config:
    #     args = read_yaml(args)
    #     pargs = sanitize_visual(args)
    # else:
    #     parse = argparse.ArgumentParser()
    #     args = parse.parse_args()
    #     args.config = [recipename]
    #     args = read_yaml(args)
    #     pargs = False
    #
    # if pargs:
    #     templatenames = wildcard(pargs.templatename)
    #     cplxnames = wildcard(pargs.cplxname)
    #     recipe = generate_recipe(pargs.fms)
    # else:

    templatenames = wildcard(templatename)
    cplxnames = wildcard(cplxname)
        # recipe = generate_recipe(args)

    all_fms = {}

    for _fm_ in recipe:
        all_distances = []
        fm_distances = []

        for ff, file in enumerate(cplxnames):

            # get cog of template and center about (0,0,0)

            template = read_pdb(templatenames[ff])

            #
            # cog of every residue
            #

            template_cogs = []

            for mol in template:
                stepx = []
                stepy = []
                stepz = []
                for pt in mol['coords']:
                    stepx.append(pt[0])
                    stepy.append(pt[1])
                    stepz.append(pt[2])

                px = np.sum(stepx) / len(stepx)
                py = np.sum(stepy) / len(stepy)
                pz = np.sum(stepz) / len(stepz)

                template_cogs.append([px, py, pz])

            #
            # cog of every residue
            #

            ampsa_cogs = []
            fm_cogs = []

            mols = read_pdb(file)
            fl = []
            for _i, mm in enumerate(mols):

                if mm['name'] == enfm[_fm_]:

                    stepx = []
                    stepy = []
                    stepz = []

                    for (x, y, z) in mm['coords']:
                        stepx.append(x)
                        stepy.append(y)
                        stepz.append(z)

                    px = np.sum(stepx) / len(stepx)
                    py = np.sum(stepy) / len(stepy)
                    pz = np.sum(stepz) / len(stepz)

                    ampsa_cogs.append([px, py, pz])
                    fm_cogs.append([px, py, pz])
                else:
                    stepx = []
                    stepy = []
                    stepz = []

                    for (x, y, z) in mm['coords']:
                        stepx.append(x)
                        stepy.append(y)
                        stepz.append(z)

                    px = np.sum(stepx) / len(stepx)
                    py = np.sum(stepy) / len(stepy)
                    pz = np.sum(stepz) / len(stepz)

                    fm_cogs.append([px, py, pz])

            # flatten the full pairwise-distance matrix rather than keeping
            # only the row for the first template residue / first complex
            # molecule, which silently discarded every other one
            dists = cdist(template_cogs, ampsa_cogs)
            dists = dists.flatten().tolist()

            adists = cdist(fm_cogs, ampsa_cogs)
            adists = adists.flatten().tolist()

            all_distances.extend(dists)
            fm_distances.extend(adists)

        all_distances = np.array(all_distances) / 10.
        fm_distances = np.array(fm_distances) / 10.
        all_fms[_fm_] = [all_distances, fm_distances]

    trajectory_kde(all_fms)
    return all_fms

def trajectory_kde(distances):
    nature = LinearSegmentedColormap.from_list('nature', ["#FEEFF9", "#825CA6"], N=100)

    fig, ax = plt.subplots(1, len(distances), sharex=True, sharey=True, figsize=(3*len(distances), 4))

    sns.set_theme(style="white")

    for i, _fm_ in enumerate(distances):
        ax[i].set_xlim([0, 10])
        ax[i].set_ylim([0, 10])
        ax[i].set(aspect='equal')
        ax[i].set_xticks([0, 2, 4, 6, 8, 10])
        ax[i].set_yticks([0, 2, 4, 6, 8, 10])

        ax[i].set_title(_fm_)

        ndata = {'x': distances[_fm_][0], 'y': distances[_fm_][1]}

        ndf = pd.DataFrame(ndata)

        sns.kdeplot(ax=ax[i], x=ndf.x, y=ndf.y, data=ndf,
                    fill=True, cmap=nature)

        ax[i].set_xlabel('FM-Residue Distance (nm)')

    ax[0].set_ylabel('FM-FM Distance (nm)')

    add_headers(fig, row_headers=['No Polymerization', "Polymerization"], fontsize=14)

    fig.tight_layout()
    # fig.subplots_adjust(right=0.8)
    cbar = fig.colorbar(cm.ScalarMappable(cmap=nature), ax=ax, #fraction=0.046, pad=0.03,
                        ticks=[0, 1])  # .ravel().tolist()
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_yticklabels(["0", "1"])
    cbar.set_label('Spatial Probability')

    plt.show()