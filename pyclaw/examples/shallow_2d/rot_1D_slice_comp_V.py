#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np


times = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
# get the data to plot
fort_output1 = './SRD_examples/Linear/_output/_gauges/gauge0.5_0.8.txt'
# fort_output1 = './SRD/General angle/Obtuse V/_output/_gauges/gauge0.5_0.3.txt'
# fort_output2 = './hbox_gen_angle/_output/_gauges/gauge0.5855_0.26399.txt'#0.5855_0.26399  0.41451_0.7169
# fort_output2 = './geoclaw_res/ObtuseV/OT/gauge00002.txt'
fort_output2 = './Mapped/_output_Lmapped/_gauges/gauge0.5_0.775.txt'#41.txt'


h_1 = np.zeros((10,10))
h_2 = np.zeros((10,10))  # column is the slice, and rows show the time progression.
plt.figure()
for i in range(1):
    # num = str(i+1)
    # fort_output_1 = fort_output1 + num + ".txt"
    # fort_output_2 = fort_output2 + num + ".txt"
    output1 = open(fort_output1,'r')
    output2 = open(fort_output2,'r')
    # output3 = open(fort_output3,'r')
    # get the right time step output
    # extract the h:
    t = []
    t2 = []
    # t3 = []
    h = []
    h2 = []
    # h3 = []
    for line in output1:
        #if even%2 == 0 and k != num_cells:
        line = line.rstrip("\n")
        hhu1 = line.split()
        # if (abs(float(hhu1[0]))>1.2):
        #     break
        # for j in range(7):
            # if (abs(float(hhu1[0])-times[j])<1e-8):
        t.append(float(hhu1[0]))
        h.append(float(hhu1[1])) #1
    # h_1[i,:] = np.asarray(h)
    for line in output2:
        line = line.rstrip("\n")
        hhu2 = line.split()
        # if (abs(float(hhu2[0]))>1.2):
        #     break
        # for j in range(7):
        #     if (abs(float(hhu2[0])-times[j])<1e-8):
        t2.append(float(hhu2[0]))
        h2.append(float(hhu2[1])) #+0.28125) #2 for geoclaw
    # for line in output3:
    #     #if even%2 == 0 and k != num_cells:
    #     line = line.rstrip("\n")
    #     hhu1 = line.split()
    #     if (abs(float(hhu1[0]))>1.):
    #         break
    #     # for j in range(7):
    #         # if (abs(float(hhu1[0])-times[j])<1e-8):
    #     t3.append(float(hhu1[1]))
    #     h3.append(float(hhu1[2]))
    # plt.figure()
    # plt.ylim(1.1,1.4)
    plt.plot(t,h,'k-',label="SRD")
    plt.plot(t2,h2,'r:',label="Mapped")
    # plt.plot(t3,h3,'b-.',label='GeoClaw')


    # plt.title("1D Slice at t="+str(times[i]))
    plt.ylabel('Height')
    plt.xlabel('Time')
    # plt.legend(loc='upper right')

    # h_2[i,:] = np.asarray(h2)
# plt.annotate('Gauge '+str(i), xy=(0.59, 1.003), xytext=(0.5, 1.017),
#         arrowprops=dict(arrowstyle="->",
#                         connectionstyle="angle3,angleA=0,angleB=-90"))
# plt.annotate('Gauge '+str(i-1), xy=(0.5, 1.05), xytext=(0.5, 1.09),
#         arrowprops=dict(arrowstyle="->",
#                         connectionstyle="angle3,angleA=0,angleB=-90"))
# plt.annotate('Gauge '+str(i+1), xy=(0.6, 1.0), xytext=(0.5,0.992),
#         arrowprops=dict(arrowstyle="->",
#                         connectionstyle="angle3,angleA=0,angleB=90"))
# plt.plot(np.linspace(0,1.4,10),np.ones(10)*1.5,'k:')
plt.legend(loc='upper left')
plt.savefig("SRD_LOT0508")
# print(h1[0,:],hu1[0,:])

# for i in range(7):
#     plt.figure()
#     # plt.ylim(0,2.6)
#     plt.plot(range(1,8),h_1[:,i],'k*',label="Rotated equivalent")
#     plt.plot(range(1,8),h_2[:,i],'ro',label="Rotated")
#     plt.title("1D Slice at t="+str(times[i]))
#     plt.ylabel('Height')
#     plt.xlabel('Gauge locations')
#     plt.legend(loc='upper right')
#     plt.savefig("1dfig"+str(i+1))
