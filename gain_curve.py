#!/usr/bin/env python
# -*- coding: utf-8 -*-

import plotly.offline as py
import plotly.graph_objs as go
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('results/gain.csv')
data.head()

#************************ PLOTLY ************************

x = [1, 2, 3, 4, 5, 6]
x_rev = x[::-1]

# NCX+Ca
y1 = data['NCX+Ca'].values
y1_upper = []
y1_lower = []
for i in y1:
    y1_upper.append(i+1)
    y1_lower.append(i-1)
y1_lower = y1_lower[::-1]

# Ca
y2 = data['Ca'].values
y2_upper = []
y2_lower = []
for i in y2:
    y2_upper.append(i+1)
    y2_lower.append(i-1)
y2_lower = y2_lower[::-1]

# NCX+Ca+Na
y3 = data['NCX+Ca+Na'].values
y3_upper = []
y3_lower = []
for i in y3:
    y3_upper.append(i+1)
    y3_lower.append(i-1)
y3_lower = y3_lower[::-1]

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=x+x_rev,
    y=y1_upper+y1_lower,
    fill='toself',
    fillcolor='rgba(0,100,80,0.2)',
    line_color='rgba(255,255,255,0)',
    showlegend=False,
    name='NCX + Ca inter',
))
fig.add_trace(go.Scatter(
    x=x+x_rev,
    y=y2_upper+y2_lower,
    fill='toself',
    fillcolor='rgba(0,176,246,0.2)',
    line_color='rgba(255,255,255,0)',
    name='Ca inter',
    showlegend=False,
))
fig.add_trace(go.Scatter(
    x=x+x_rev,
    y=y3_upper+y3_lower,
    fill='toself',
    fillcolor='rgba(240,52,52,0.2)',
    line_color='rgba(255,255,255,0)',
    showlegend=False,
    name='NCX + Ca + Na',
))
fig.add_trace(go.Scatter(
    x=x, y=y1,
    line_color='rgb(0,100,80)',
    name='NCX + Ca',
))
fig.add_trace(go.Scatter(
    x=x, y=y2,
    line_color='rgb(0,176,246)',
    name='Ca inter',
))
fig.add_trace(go.Scatter(
    x=x, y=y3,
    line_color='rgb(240,52,52)',
    name='NCX + Ca + Na',
))

# Edit the layout
fig.update_layout( xaxis_title='Range (#cells)',
                   yaxis_title='Gain (dB)',
                   legend=dict(x=.68, y=0.45))


fig.update_traces(mode='lines+markers')
fig.show()

# fig = go.Figure(data=data, layout=layout)
# py.plot(fig, filename='delay_chart.html')

#************************ PLOTLY ************************


# data = pd.read_csv('erros.csv')
# data.head()

# E1 = data['ALPHA_MLP'].values
# E2 = data['ETA_MLP'].values
# E3 = data['ETA_ALPHA_MLP'].values
# E4 = data['NORMAL_MLP'].values
# E5 = data['MLP_ONLINE'].values
# E6 = data['RPROP'].values
# E7 = data['RBF_BP'].values
# E8 = data['RBF_BATCH'].values
# X = np.array(list(zip(E1/max(E1), E2/max(E2), E3/max(E3), E4/max(E4), E5/max(E5), E6/max(E6), E7, E8)), dtype=np.float32)


# epocas = np.linspace(0, X.shape[0], X.shape[0])
# legendas = ['MLP com Alpha', 'MLP com Eta variável', 'MLP com Eta variável e Alpha', 'MLP Original', 'MLP Online', 'RPROP', 'RBF Original', 'RBF Batch']
# estilos_linha = ['-', '-', '-', '--', '-.', '-', '-', '-']

# fig, ax = plt.subplots()
# """
# for i in range(X.shape[-1]):
#     ax.semilogy(epocas, X[:, i], label=legendas[i], linestyle=estilos_linha[i])
# """
# ax.semilogy(epocas, X[:, 6])

# ax.set(xlabel='Épocas de Treinamento', ylabel='Eav')
# ax.grid()
# ax.legend(loc=1)
# fig.savefig("Erros.png")
# plt.show()
