import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DIM_X = 5 #Pegar com base no header
DIM_Y = 3


data = pd.read_csv("data.txt")

#print(data.iloc[selected_line:selected_line+1,:].values[0])


for selected_line in range(0, data.shape[0]):

    npdata = np.array([])
    for i in range(0, DIM_X):
        for j in range(0, DIM_Y):
            npdata = np.append(npdata,data.iloc[selected_line:selected_line+1,DIM_Y*i+j].values[0])

    npdata = npdata.reshape(DIM_X,DIM_Y)
    print(npdata)


    if selected_line == 0:
        plt.title("Diffusion")
        p = plt.imshow(npdata)
        plt.colorbar(p)
        plt.show()
        p = plt.imshow(npdata)
        plt.colorbar(p)
        fig = plt.gcf()
        plt.clim()   # clamp the color limits
        plt.title("Diffusion")
    else:
        p.set_data(npdata)

    print("step", selected_line)
    plt.pause(0.1)
plt.show()
