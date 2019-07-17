import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#DIM_X = 5 #Pegar com base no header
#DIM_Y = 3


data = pd.read_csv("Temp/data.txt")

#Finding the tissue's dimensions
underline = 0;
a = data.columns[-1]
aux_X = ''
aux_Y = ''
for i in a:
    if i == '_':
        underline = 1
    else:
        if underline == 0:
            aux_X += i;
        else:
            aux_Y += i
DIM_X = int(aux_X) + 1
DIM_Y = int(aux_Y) + 1

#print(data.iloc[selected_line:selected_line+1,:].values[0])


for selected_line in range(0, data.shape[0]):

    npdata = np.array([])
    for i in range(0, DIM_X):
        for j in range(0, DIM_Y):
            npdata = np.append(npdata,data.iloc[selected_line:selected_line+1,DIM_Y*i+j].values[0])

    npdata = (npdata.reshape(DIM_X,DIM_Y)).transpose()
    print(npdata)


    if selected_line == 0:
        plt.title("Tissue - Initial State")
        p = plt.imshow(npdata, cmap = 'hot', vmax = 6, vmin = 0.1)
        plt.colorbar(p)
        print("\nInitial state of the tissue. Close the plot window to proceed.\n")
        plt.show()
        p = plt.imshow(npdata, cmap = 'hot', vmax = 6, vmin = 0.1)
        plt.colorbar(p)
        fig = plt.gcf()
        plt.clim()   # clamp the color limits
        plt.title("Reactions on the Tissue")
    else:
        p.set_data(npdata)

    print("\nStep", selected_line)
    plt.pause(0.1)
print("\nThe simulation ended. Close the plot window to finish the application.\n\n")
plt.show()
