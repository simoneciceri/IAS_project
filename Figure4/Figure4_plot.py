import matplotlib.pyplot as plt 
import matplotlib
import seaborn as sns
import numpy as np
import math

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


plot_colorbar = False
bool_isolines = True



dir_list = [    
    "globalfeedback",
    "w12",
    "w13",
    "w23",
    "isolatePPC",
    "isolatePFC"
    ]

titlelist = [   
    "Global feedback",
    r"PPC $\to$ V1",
    r"PFC $\to$ V1",
    r"PFC $\to$ PPC",
    "Isolate PPC",  #isolate2
    "Isolate PFC",  #isolate3
    ]


x_lims = {                              #not used anymore
    "global_feedback":[0.6, 1.51],
    "w21":[0.4, 1.51],
    "w31":[0., 1.51],
    "w32":[0.6, 1.51],
    "isolate2":[0.9, 1.1],
    "isolate3":[0.8, 1.3]
    }

xticks=[0., 0.5, 1., 1.5]
yticks=[0., 1, 2, 3, 4]



NA = 61                                         # number of homotopy param values 
NI = 100                                        # number of applied current values
NIC = 50                                        # number of different initial conditions for each value of current and alpha
alpha_list = np.arange(0., 1.51, 0.025)         # homotopy parameter values
maxcurrent = 4
currents = np.arange(0, maxcurrent, 0.04)       # applied current values

title_idx = 0

# FOR EACH ENSEMBLE OF CONNECTIONS
for directory in dir_list:
    prob_fb = np.zeros([NA, NI])                   # probability of feedback bumps
    prob_sb = np.zeros([NA, NI])                   # probability of only one bump
    prob_ov = np.zeros([NA, NI])                   # probability of overshooting

# FOR EACH VALUE OF HOMOTOPY PARAMETER ALPHA
    for ai in range(NA):

        alpha = alpha_list[ai]                                              #alpha            
        filename = "./" + directory + "/v1Int_a{:g}.dat".format(alpha)      #open file
        f = open(filename, "r")
        lines = f.readlines()
        f.close()

        #reading & assembling the distribution of integrals
        distr_int = []
        for line in lines:
            distr_int += [np.array( eval(line)) ]
        distr_int = np.stack(distr_int)
        #print(distr_int.shape)

        #findind the probability of second bump for each value of applied current
        for ni in range(NI):
            hit=0
            for nic in range(NIC):
                if(distr_int[ni, nic]<20): prob_sb[ai, ni]+=1
                elif(distr_int[ni, nic]>20 and distr_int[ni, nic]<35): prob_fb[ai,ni]+=1
                else: prob_ov[ai, ni]+=1
            prob_fb[ai, ni] /= NIC
            prob_sb[ai, ni] /= NIC
            prob_ov[ai, ni] /= NIC

    prob_fb = np.transpose(prob_fb)
    prob_sb = np.transpose(prob_sb)
    prob_ov = np.transpose(prob_ov)


    # ========= HEATMAP OF THE FEEDBACK BUMP PROBABILITY ===========
    plt.contourf(alpha_list, currents, prob_fb, cmap="viridis", levels=50, vmin=0.01, vmax=1, 
                        norm=matplotlib.colors.PowerNorm(0.4))
                        #normalizing the colorbar with power Law weight
    
    if plot_colorbar:
        cbar = plt.colorbar(ticks=[0,0.25, 0.5, 0.75, 1])
        cbar.set_label('Feedback-bump Probability',size=20)
        cbar.ax.tick_params(labelsize=15)

     
    if(bool_isolines):
        # ========= ISOMETRYC LINES OF THE OVERSHOOTING PROB: 99% ===========
        plt.contour(alpha_list, currents, prob_ov, levels=[0.99], colors="tab:orange", linestyles="-")

        # ========= ISOMETRYC LINES OF THE SINGLE-BUMP PROB: 99% ===========
        plt.contour(alpha_list, currents, prob_sb, levels=[0.99], colors="tab:green", linestyles="-")


    # some cosmetics
    plt.rcParams["hatch.linewidth"] = 4
    plt.xlim([0, 1.5])
    plt.ylim(0, 3.9)
    plt.xlabel("Morphing parameter", fontsize=20)
    plt.ylabel("Maximum applied current (mA)", fontsize=19)
    #plt.ylabel(r"Maximum applied current $I_{max}$ (mA)", fontsize=19)
    #plt.title("{}".format(titlelist[title_idx]), fontsize=22)
    plt.xticks(ticks=xticks, fontsize=18)
    plt.yticks(ticks=yticks, fontsize=18)
    plt.text(0.1, 0.5, "{}".format(titlelist[title_idx]), fontsize=22, color="white")

    #change figure size
    fig = plt.gcf()
    [figWidth,figHeight] = plt.rcParams.get('figure.figsize')
    fig.set_figheight(figHeight + 0.7)

    # save figures in SVG format
    fig.savefig("./subfigures/{}.pdf".format(dir_list[title_idx]))
    title_idx +=1

    plt.show()
    #break

