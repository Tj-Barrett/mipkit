import matplotlib.pyplot as plt
import numpy as np
import os

import matplotlib.animation as animation


def polyanimation(args):
    def update(frame):
        x = time[frame : frame + 1]
        y = mwnavg[frame : frame + 1]
        y1 = mwwavg[frame : frame + 1]
        y2 = npoly[frame : frame + 1]
        pt.set_xdata(x)
        pt.set_ydata(y)
        pt1.set_xdata(x)
        pt1.set_ydata(y1)
        pt2.set_xdata(x)
        pt2.set_ydata(y2)

    basedir = args.basedir[0] if isinstance(args.basedir, list) else args.basedir

    if args.xstats:
        xstats = os.path.join(basedir, args.xstats[0])
    else:
        print(
            "No polyermization stats file specified. Defaulting to work/stats_polymerization ..."
        )
        xstats = os.path.join(basedir, "work/stats_polymerization.txt")

    if args.mwstats:
        mwstats = os.path.join(basedir, args.mwstats[0])
    else:
        print(
            "No molecular weight stats file specified. Defaulting to work/stats_mw_distributions ..."
        )
        mwstats = os.path.join(basedir, "work/stats_mw_distributions.txt")

    if float(args.giftime[0]):
        giftime = args.giftime[0] * 1000.0  # gif time in ms
    else:
        print("Gif time is not a number. Quiting ...")
        quit()

    if args.gifname:
        gifname = os.path.join(basedir, args.gifname[0])
    else:
        print("No gifname given. Defaulting.")
        gifname = os.path.join(basedir, "default.gif")

    step = []
    nmono = []
    npoly = []
    mwnavg = []
    mwwavg = []
    mwmax = []
    mwmin = []

    time = []

    totaltime = giftime  # gif time in ms

    with open(xstats, "r") as file:
        lines = file.readlines()
    for _l, l in enumerate(lines):
        if _l < 3:
            continue
        else:
            sl = l.split()
            # print(sl)
            step.append(float(sl[0]))
            time.append(0.5 * float(sl[0]))
            nmono.append(float(sl[1]))
            npoly.append(float(sl[2]))
            mwnavg.append(float(sl[3]))
            mwwavg.append(float(sl[4]))
            mwmax.append(float(sl[5]))
            mwmin.append(float(sl[6]))

    mws = []
    with open(mwstats, "r") as file:
        lines = file.readlines()
    for _l, l in enumerate(lines):
        if _l < 3:
            continue
        else:
            sl = l.split()
            ns = []
            for n in sl[1:]:
                if float(n) > 0:
                    ns.append(float(n))
            mws.append(ns)

    fig, ax = plt.subplots(nrows=1, ncols=1)

    ax1 = ax  # [0]

    ax1.plot(time, mwnavg, "b")
    ax1.plot(time, mwwavg, "b:")
    pt = ax1.plot(time[0], mwnavg[0], "bo")[0]
    pt1 = ax1.plot(time[0], mwwavg[0], "bo")[0]
    ax1.set_xlabel("Time [ns]")
    ax1.set_ylabel("Molecular Weight [g/mol]")
    ax1.set_ylim([0, int(np.ceil(np.max(mwwavg) / 100.0) * 100)])
    ax1.tick_params(axis="y", colors="blue")
    ax1.spines["left"].set_color("blue")
    ax1.set_xlim([0, 25])

    ax2 = ax1.twinx()
    ax2.plot(time, npoly, "r")
    pt2 = ax2.plot(time[0], npoly[0], "ro")[0]
    ax2.set_ylabel("Polymer Count")
    ax2.set_ylim([0, int(np.ceil(np.max(npoly) / 10.0) * 10)])
    ax2.tick_params(axis="y", colors="red")
    ax2.spines["right"].set_color("red")

    ax1.legend(["Number Average", "Weight Average"])  # , 'Polymer Count')

    interval = totaltime / len(step)

    ani = animation.FuncAnimation(
        fig=fig, func=update, frames=len(step), interval=interval
    )
    plt.show()
    ani.save(gifname)
