import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DIM_X = 5 #Pegar com base no header
DIM_Y = 3

selected_line = 188 #depois posso criar um iterador

data = pd.read_csv("data.txt")

#print(data.iloc[selected_line:selected_line+1,:].values[0])

npdata = np.array([])

for i in range(0, DIM_X):
    for j in range(0, DIM_Y):
        npdata = np.append(npdata,data.iloc[selected_line:selected_line+1,DIM_Y*i+j].values[0])

npdata = npdata.reshape(DIM_X,DIM_Y)
print(npdata)


y = np.arange(0,DIM_Y)
plt.xticks(y)
plot = plt.imshow(npdata, cmap='viridis', vmin=0,vmax=0.5)
plt.colorbar(plot)
plt.tight_layout()
plt.show()

#plotar uma animacao do heatmap se atualizando
