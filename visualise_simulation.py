#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Oct 24 16:35:49 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

- This code can be used to visualise the synthetic vascular network, along with 
the haematocrit in each vessel and the resultant oxygen concentration.

- Set only_hypoxia=1 if you just want to visualise the hypoxia (red is 
normoxic, blue is hypoxic) or only_hypoxia=0 to see the oxygen field itself.

- Set mark_nodes=1 if you want to see the nodes or mark_nodes=0 to just see the 
edges.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Import custom functions
from convert_oxygen_units import *
from scipy.ndimage import zoom

# Set LaTex-style font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22})
import matplotlib as mpl
mpl.rcParams['axes.grid'] = False  # Disable grid globally

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to visualise a simulation
def visualise_simulation(G_filtered, only_hypoxia=0, mark_nodes=0, hypoxic_threshold=4500):

    # Extract the required data
    O2_field = G_filtered.graph['oxygen_field']
    NF = G_filtered.graph['NF'] 

    # Set the simulation boundaries for plotting
    x_min = 0
    x_max = 1105.0
    y_min = 0
    y_max = 1105.0 + (2.0*85.0)

    # Create a meshgrid for the new dimensions based on the segment boundaries
    new_x_size = int(x_max - x_min) + 1  # Width of the new image
    new_y_size = int(y_max - y_min) + 1  # Height of the new image

    # Stretch the O2_field to match the new dimensions
    zoom_x = new_x_size / O2_field.shape[1]
    zoom_y = new_y_size / O2_field.shape[0]
    O2_field_stretched = zoom(O2_field, (zoom_y, zoom_x))
     
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.grid(False)
        
    # Determine min and max haematocrit values across all edges
    haematocrit_values = [edge_data.get('haematocrit', 0) for _, _, edge_data in G_filtered.edges(data=True)]
    haematocrit_min = 0
    haematocrit_max = max(haematocrit_values)
    
    # Overlay the segments extracted from G_filtered, colouring by haematocrit
    for edge in G_filtered.edges(data=True):
        # Extract edge information
        edge_data = edge[2]
        if 'segment_coordinates' in edge_data:
            x1, y1, x2, y2 = edge_data['segment_coordinates']
            haematocrit = edge_data.get('haematocrit', 0)  # Default to 0 if haematocrit not found
            color = plt.cm.jet((haematocrit - haematocrit_min) / (haematocrit_max - haematocrit_min))
            ax.plot([x1, x2], [y1, y2], c=color, lw=3)

    # Add a color bar to indicate haematocrit levels, with range based on min and max haematocrit values
    sm = plt.cm.ScalarMappable(cmap='jet', norm=plt.Normalize(vmin=haematocrit_min, vmax=haematocrit_max))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.02, fraction=0.05, location='right')
    cbar.ax.set_position([cbar.ax.get_position().x0, cbar.ax.get_position().y0, cbar.ax.get_position().width, cbar.ax.get_position().height * 0.45])
    cbar.set_label('haematocrit (concentration)')

    # Add markers for nodes
    if mark_nodes == 1:
        for node in G_filtered.nodes():
            connected_edges = list(G_filtered.edges(node, data=True))
            if connected_edges:
                
                # Use the coordinates from the first connected edge involving the node
                edge_data = connected_edges[0][2]
                if 'segment_coordinates' in edge_data:
                    x1, y1, x2, y2 = edge_data['segment_coordinates']
                    
                    # Determine which endpoint is the current node
                    if node == connected_edges[0][0]:
                        x, y = x1, y1
                    else:
                        x, y = x2, y2
                    ax.scatter(x, y, c='black', s=50, zorder=3)  # Add a black marker for each node

    # Display the stretched O2 field with custom colouring
    if only_hypoxia == 1:
        O2_field_coloured = np.where(O2_field_stretched < convert_nM_to_mmHg(hypoxic_threshold), 0, 1) 
        im = ax.imshow(O2_field_coloured, cmap='coolwarm', origin='lower', extent=[x_min, x_max, y_min, y_max], vmin=0, vmax=1, alpha=0.5)
    elif only_hypoxia == 0:
        O2_field_coloured = O2_field_stretched
        im = ax.imshow(O2_field_stretched, cmap='coolwarm', origin='lower')
        
        # Add a second color bar for O2 field if only_hypoxia is 0
        sm_o2 = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=O2_field_stretched.min(), vmax=O2_field_stretched.max()))
        sm_o2.set_array([])
        cbar_o2 = fig.colorbar(sm_o2, ax=ax, orientation='vertical', pad=0.02, fraction=0.05, location='right')
        cbar_o2.ax.set_position([cbar.ax.get_position().x0, cbar_o2.ax.get_position().y0 + 0.41, cbar_o2.ax.get_position().width, cbar_o2.ax.get_position().height * 0.45])
        cbar_o2.set_label('oxygen (mmHg)')

    # Set labels and title
    ax.set_xlabel('x-axis (μm)')
    ax.set_ylabel('y-axis (μm)')
    ax.set_title(f'NF = {NF:.2f}')  # Use the NF value for title
    plt.show()

# =============================================================================
# MAIN 
# =============================================================================

# Enter the details of the simulation data to analyse
which_selection = 33  # Ranges from 1 to 100
n_kills = 39

# Load the graph            
load_path = f'data/network_s{which_selection}_k{n_kills}.pkl'
with open(load_path, 'rb') as f:
    loaded_network = pickle.load(f)

# Visualise the graph and the resultant oxygenation
visualise_simulation(loaded_network, only_hypoxia=1, mark_nodes=1)
    