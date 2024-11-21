# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 20:28:11 2024

@author: domin
"""

import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go

with open("best.txt", 'r') as file:
    file.readline()
    max_compactness = file.readline()
    _, _, max_compactness, _, _, _ = max_compactness.split()
max_compactness = float(max_compactness)

def read_file(path):
    with open(path, 'r') as file:
        energy = file.readline()
        _, energy = energy.split()
        energy = float(energy)
        Compactness = file.readline()
        _, Compactness = Compactness.split()
        Compactness = float(Compactness) / max_compactness
        temperature = file.readline()
        _, temperature = temperature.split()
        temperature = float(temperature)
    
    df = pd.read_csv(path, sep=r'\s+', names=["L", "x", "y", "z"], skiprows=3)
    letters = df["L"].to_numpy()
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()
    z = df["z"].to_numpy()
    
    return letters, x, y, z, energy, Compactness

def map_amino_acid_to_colors(amino_acids):
    amino_acid_colors = {
        'A': 'blue', 'V': 'cyan', 'L': 'skyblue', 'I': 'deepskyblue', 
        'M': 'dodgerblue', 'F': 'royalblue', 'W': 'steelblue', 
        'Y': 'lightblue', 'P': 'lightcyan', 'C': 'lightseagreen', 
        'G': 'lightsteelblue', 'S': 'red', 'T': 'orange', 
        'N': 'yellow', 'Q': 'gold', 'D': 'darkorange', 
        'E': 'darkred', 'K': 'darkred', 'R': 'firebrick', 
        'H': 'darkorange'
    }
    
    return [amino_acid_colors[aa] for aa in amino_acids], amino_acid_colors

folder_path = 'conformations'
files = sorted(os.listdir(folder_path), key=lambda x: int(x.split('.')[0]))

fig = go.Figure()
letters, x, y, z, energy, compactness = read_file(os.path.join(folder_path, files[0]))
color_values, amino_acid_colors = map_amino_acid_to_colors(letters)

# Add initial 3D trace
fig.add_trace(go.Scatter3d(
    x=x, y=y, z=z, mode='markers+lines',
    marker=dict(size=4, color=color_values),
    line=dict(color='black', width=0.5),
    name='Protein Structure'
))

# Add legend items for amino acid colors
for aa, color in amino_acid_colors.items():
    fig.add_trace(go.Scatter3d(
        x=[None], y=[None], z=[None],
        mode='markers',
        marker=dict(symbol='circle', color=color, size=10),
        name=f'{aa}',
        showlegend=True
    ))


fixed_camera = dict(
    eye=dict(x=1.25, y=1.25, z=1.25),  # Adjust x, y, z as needed
    center=dict(x=0, y=0, z=0),
    up=dict(x=0, y=0, z=1)
)

# Create frames for animation
frames = []
for i, file_name in enumerate(files):
    letters, x, y, z, energy, compactness = read_file(os.path.join(folder_path, file_name))
    color_values, _ = map_amino_acid_to_colors(letters)
    frames.append(go.Frame(
        data=[go.Scatter3d(
            x=x, y=y, z=z, mode='markers+lines',
            marker=dict(size=4, color=color_values),
            line=dict(color='black', width=0.5)
        )],
        name=str(i),
        layout=go.Layout(
            annotations=[dict(
                x=0.02, y=0.98,
                xref="paper", yref="paper",
                text=f"Energy: {energy}, compactness: {compactness:.2f}",
                showarrow=False,
                font=dict(size=14, color="black"),
                align="left"
            )]
        )
    ))

fig.frames = frames

# Add play button and frame-navigation slider
fig.update_layout(
    scene=dict(
        xaxis=dict(title="X", range=[-100, 100]), 
       yaxis=dict(title="Y", range=[-100, 100]), 
       zaxis=dict(title="Z", range=[-100, 100]),
       camera=fixed_camera,
       aspectmode="manual", 
       aspectratio=dict(x=1, y=1, z=1) 
    ),
    title="Animated 3D Scatter Plot",
    showlegend=True,
    updatemenus=[dict(
        type="buttons",
        showactive=True,
        buttons=[
            dict(
                label="Play",
                method="animate",
                args=[None, dict(frame=dict(duration=200, redraw=True), fromcurrent=True)]
            ),
            dict(
                label="Pause",
                method="animate",
                args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate", fromcurrent=True)]
            )
        ]
    )],
    sliders=[{
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {"prefix": "Frame: ", "font": {"size": 20}, "visible": True},
        "steps": [
            {
                "args": [[str(i)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                "label": str(i),
                "method": "animate"
            }
            for i in range(len(files))
        ]
    }]
)
fig.show()

fig.write_html("compactness3D2.html")