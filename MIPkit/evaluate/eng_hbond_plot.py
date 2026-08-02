import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np

from MIPkit.utils.utils import (
    moving_average,
)

def plot_eng_hbond(args):

    font1 = {'family':'serif','color':'black','size':20}
    font2 = {'family':'serif','size':20}
    fig, axs = plt.subplots(2,1, height_ratios=[3,1], figsize=(10, 8), dpi=80)

    colors = ['r','b']
    _i=-1

    mmip = "BRCA1"

    for _ in [0,1]:
        total_eng = []
        total_hbond = []
        _i+=1

        for j in [1,2,3,4,5]:
            try:
                if _i == 0:
                    ff = f'MIP-Interactions/{mmip}-MD-{j}-MIP'
                else:
                    ff = f'MIP-Interactions/{mmip}-ML-{j}-MIP'

                #
                # Interaction Energy
                #

                toteng = 0.

                for _j, f in enumerate([ff+'/epitope_lj.xvg',ff+'/epitope_coul.xvg' ]):
                    engav = []
                    _f = open(f)
                    lines = _f.readlines()[24:-1]
                    _f.close()

                    time = [0]*len(lines)
                    eng  = [0]*len(lines)

                    if _i == 0:
                        mipeng = np.array(eng)
                        nipeng = np.array(eng)
                    if _j == 0:
                        neteng = np.zeros(len(eng))

                    for i, _l in enumerate(lines):
                        _sl = _l.split()
                        time[i] = float(_sl[0])/1000.
                        eng[i] = float(_sl[1])
                        if float(_sl[0])/1000. > 30:
                            engav.append(float(_sl[1]))

                    neteng+=np.array(eng)
                    toteng+=np.mean(engav)
                total_eng.append(toteng)
                print(f'{ff} Mean : {round(toteng,3)} kJ/mol')

                # MA
                _eng = moving_average(neteng,20)
                _time = moving_average(time,20)
                #
                #
                #
                if _i == 0:
                    md, = axs[0].plot(-_time, -_eng, colors[_i], alpha=1)
                elif _i == 1:
                    ml, = axs[0].plot(-_time, -_eng, colors[_i], alpha=1)
                # Full
                axs[0].plot(_time, _eng, colors[_i], alpha=0.4)
                axs[0].plot(time, neteng, colors[_i], alpha=0.1)

                #
                # H Bonds
                #
                f = ff+"/hbond_nums.xvg"
                _f = open(f)
                lines = _f.readlines()[25:-1]
                _f.close()

                time = [0]*len(lines)
                count  = [0]*len(lines)
                cutcount = [0]*len(lines)

                avghbond = []

                for i, _l in enumerate(lines):
                    _sl = _l.split()
                    time[i] = float(_sl[0])/1000.
                    count[i] = float(_sl[1])
                    cutcount[i] = float(_sl[2])
                    if float(_sl[0])/1000. > 30:
                        avghbond.append(float(_sl[2]))

                _count = moving_average(cutcount,20)
                _time = moving_average(time,20)

                axs[1].plot(_time, _count, colors[_i], alpha=0.4)
                axs[1].plot(time, cutcount, colors[_i], alpha=0.1)

                total_hbond.append(np.mean(avghbond))

                print(f'{ff} Mean : {round(np.mean(avghbond),3)} H-Bonds')

            except (ValueError, TypeError, OSError) as e:
                print(f'Warning: skipping {ff} due to {type(e).__name__}: {e}')
                continue

        try:
            print(f'\nOverall Mean : {round(np.mean(total_eng),3)} kJ/mol')
            print(f'Overall STD : {round(np.std(total_eng),3)} kJ/mol\n')
        except (ValueError, TypeError) as e:
            print(f'Warning: could not compute overall energy statistics: {e}')
            continue
        try:
            print(f'Overall Mean : {round(np.mean(total_hbond),3)} H-Bonds\n')
        except (ValueError, TypeError) as e:
            print(f'Warning: could not compute overall H-bond statistics: {e}')
            continue

    axs[0].set_xlim([-1, 51])
    axs[1].set_xlim([-1, 51])

    axs[0].set_ylim([-1500, 50])
    axs[1].set_ylim([-1, 25])

    axs[1].set_xlabel('Time [ns]', fontdict=font1)
    axs[0].set_ylabel('Energy [kJ/mol]', fontdict=font1)
    axs[1].set_ylabel('H-Bond Count', fontdict=font1)



    axs[0].set_xticklabels(axs[0].get_xticks(), fontdict=font1)
    axs[1].set_xticklabels(axs[1].get_xticks(), fontdict=font1)
    axs[0].set_yticklabels(axs[0].get_yticks(), fontdict=font1)
    axs[1].set_yticklabels(axs[1].get_yticks(), fontdict=font1)
    axs[0].xaxis.set_tick_params(labelsize=15)
    axs[1].xaxis.set_tick_params(labelsize=15)
    axs[0].yaxis.set_tick_params(labelsize=15)
    axs[1].yaxis.set_tick_params(labelsize=15)

    axs[0].xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    axs[0].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    axs[1].xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    axs[1].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    try:
        axs[0].legend([md, ml], ['MD', 'ML'], prop=font2)
    except (ValueError, TypeError, NameError) as e:
        print("")

    fig.align_ylabels(axs)
    plt.show()
