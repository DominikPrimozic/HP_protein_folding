# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 10:49:32 2024

@author: domin
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 10:37:11 2024

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
max_compactness=float(max_compactness)

def read_file(path):
    with open(path, 'r') as file:
        energy = file.readline()
        _, energy = energy.split()
        energy = float(energy)
        Compactness = file.readline()
        _, Compactness = Compactness.split()
        Compactness = float(Compactness)/max_compactness #this is from simulation
        temperature = file.readline()
        _, temperature = temperature.split()
        temperature = float(temperature)
    
    df = pd.read_csv(path, sep=r'\s+', names=["L", "x", "y"], skiprows=3)
    letters = df["L"].to_numpy()
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()
    
    return letters, x, y, energy, Compactness

def map_amino_acid_to_colors(amino_acids):
    amino_acid_colors = {
        # Nonpolar (cold colors)
        'A': 'blue',     # Alanine
        'V': 'cyan',     # Valine
        'L': 'skyblue',  # Leucine
        'I': 'deepskyblue', # Isoleucine
        'M': 'dodgerblue', # Methionine
        'F': 'royalblue',  # Phenylalanine
        'W': 'steelblue',  # Tryptophan
        'Y': 'lightblue', # Tyrosine
        'P': 'lightcyan', # Proline
        'C': 'lightseagreen', # Cysteine (polar with some nonpolar characteristics)
        'G': 'lightsteelblue', # Glycine (nonpolar, neutral)

        # Polar (warm colors)
        'S': 'red', # Serine
        'T': 'orange',   # Threonine
        'N': 'yellow',   # Asparagine
        'Q': 'gold', # Glutamine

        # Charged, hydrophilic (warm colors)
        'D': 'darkorange',     # Aspartic Acid
        'E': 'darkred', # Glutamic Acid
        'K': 'darkred',     # Lysine
        'R': 'firebrick', # Arginine
        'H': 'darkorange'  # Histidine
    }
    
    return [amino_acid_colors[aa] for aa in amino_acids],amino_acid_colors

folder_path = 'conformations'
files = sorted(os.listdir(folder_path), key=lambda x: int(x.split('.')[0]))

fig = go.Figure()
letters, x, y, energy, compactness = read_file(os.path.join(folder_path, files[0]))
color_values, amino_acid_colors = map_amino_acid_to_colors(letters)

# Add initial trace
fig.add_trace(go.Scatter(
    x=x, y=y, mode='markers+lines',
    marker=dict(size=4, color=color_values),
    line=dict(color='black', width=0.5),
    name='Protein Structure'
))

# Add legend items for amino acid colors
for aa, color in amino_acid_colors.items():
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', color=color, size=10),
        name=f'{aa}',
        showlegend=True
    ))

# Create frames for animation
frames = []
for i, file_name in enumerate(files):
    letters, x, y, energy, compactness = read_file(os.path.join(folder_path, file_name))
    color_values, _ = map_amino_acid_to_colors(letters)
    frames.append(go.Frame(
        data=[go.Scatter(
            x=x, y=y, mode='markers+lines',
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
    title="Animated Scatter Plot",
    xaxis=dict(title="X"),
    yaxis=dict(title="Y"),
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

fig.write_html("compactness3.html")