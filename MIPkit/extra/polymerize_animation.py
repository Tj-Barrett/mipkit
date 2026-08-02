import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation


def update(frame):
    x = time[frame : frame + 1]
    y = mwnavg[frame : frame + 1]
    y2 = npoly[frame : frame + 1]
    pt.set_xdata(x)
    pt.set_ydata(y)
    pt2.set_xdata(x)
    pt2.set_ydata(y2)


filename = (
    "/home/tj/Documents/Polymerization/cg-vina/3.5 A Template/stats_polymerization.txt"
)

step = []
nmono = []
npoly = []
mwnavg = []
mwwavg = []
mwmax = []
mwmin = []

time = []

totaltime = 5000  # gif time in ms

with open(filename, "r") as file:
    lines = file.readlines()
for _l, l in enumerate(lines):
    if _l == 0:
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


filename = "/home/tj/Documents/Polymerization/cg-vina/3.5 A Template/stats_mw_distributions.txt"

mws = []
with open(filename, "r") as file:
    lines = file.readlines()
for _l, l in enumerate(lines):
    if _l == 0:
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
pt = ax1.plot(time[0], mwnavg[0], "bo")[0]
ax1.set_xlabel("Time [ns]")
ax1.set_ylabel("Number Average MW [g/mol]")
ax1.set_ylim([0, 500])
ax1.tick_params(axis="y", colors="blue")
ax1.spines["left"].set_color("blue")
ax1.set_xlim([0, 25])

ax2 = ax1.twinx()
ax2.plot(time, npoly, "r")
pt2 = ax2.plot(time[0], npoly[0], "ro")[0]
ax2.set_ylabel("Polymer Count")
ax2.set_ylim([0, 20])
ax2.tick_params(axis="y", colors="red")
ax2.spines["right"].set_color("red")


# counts, bins = np.histogram(mws[0])
# hax[1].hist(counts, bins, color='b')
# # pt = ax1.plot(time[0], mwnavg[0], 'bo')[0]
# ax[1].set_ylabel('Count')
# ax[1].set_xlabel('Number Average MW [g/mol]')
# ax[1].set_xlim([0,500])
# ax[1].set_ylim([0,500])


interval = totaltime / len(step)

ani = animation.FuncAnimation(fig=fig, func=update, frames=len(step), interval=interval)
plt.show()
ani.save("images/35 Template.gif")


# https://stackoverflow.com/questions/51517685/combine-several-gif-horizontally-python
# import imageio

# gif1=imageio.get_reader('images/35A Teplate.gif')
# gif2=imageio.get_reader('images/mw-polymers-animation.gif')

# #If they don't have the same number of frame take the shorter
# number_of_frames = min(gif1.get_length(), gif2.get_length())

# #Create writer object
# new_gif = imageio.get_writer('output.gif')

# for frame_number in range(number_of_frames):
#     img1 = gif1.get_next_data()
#     img2 = gif2.get_next_data()
#     #here is the magic
#     new_image = np.hstack((img1, img2))
#     new_gif.append_data(new_image)

# gif1.close()
# gif2.close()
# new_gif.close()
