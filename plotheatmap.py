import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#DIM_X = 5 #Pegar com base no header
#DIM_Y = 3


dataCa = pd.read_csv("temp/data.txt")

#Finding the tissue's dimensions
underline = 0
a = dataCa.columns[-1]
aux_X = ''
aux_Y = ''
for i in a:
    if i == '_':
        underline = 1
    else:
        if underline == 0:
            aux_X += i
        else:
            aux_Y += i
DIM_X = int(aux_X) + 1
DIM_Y = int(aux_Y) + 1

#print(dataCa.iloc[selected_line:selected_line+1,:].values[0])


for selected_line in range(0, dataCa.shape[0]):

    npdataCa = np.array([])
    for i in range(0, DIM_X):
        for j in range(0, DIM_Y):
            npdataCa = np.append(npdataCa,dataCa.iloc[selected_line:selected_line+1,DIM_Y*i+j].values[0])

    npdataCa = (npdataCa.reshape(DIM_X,DIM_Y)).transpose()
    print("\n====================\n\n",npdataCa)


    if selected_line == 0:
        plt.title("Initial State - Sodium:")
        p = plt.imshow(npdataCa, cmap = 'hot', vmax = 20000, vmin = 15000)
        plt.colorbar(p)
        plt.xticks(np.arange(0, DIM_X))
        plt.yticks(np.arange(0, DIM_Y))

        print("\nTime", selected_line)
        print("\nInitial state of the tissue. Close the plot window to proceed.\n")
        plt.show()
        p = plt.imshow(npdataCa, cmap = 'hot', vmax = 20000, vmin = 15000)
        plt.colorbar(p)
        plt.xticks(np.arange(0, DIM_X))
        plt.yticks(np.arange(0, DIM_Y))
        fig = plt.gcf()
        plt.clim()   # clamp the color limits
        plt.title("Sodium diffusion on the tissue")
    else:
        p.set_data(npdataCa)
        print("\nTime", selected_line)
    plt.pause(0.1)
print("\nThe simulation ended. Close the plot window to finish the application.\n\n")
plt.show()
