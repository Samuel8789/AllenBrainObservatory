# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 12:41:05 2020

@author: sp3660
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button

raw_data = h5py.File('G:\CodeTempRawData\EnsembleMethodsAllenExc\ophys_experiment_528402271.h5', 'r')
raw_data['data']

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2, top=0.9)

ax_T = fig.add_axes([0.2, 0.95, 0.5, 0.075])
ax_T.spines['top'].set_visible(True)
ax_T.spines['right'].set_visible(True)

ax_play=fig.add_axes([0.2, 0.05, 0.5, 0.075])
ax_play.spines['top'].set_visible(True)
ax_play.spines['right'].set_visible(True)

s_T = Slider(ax=ax_T, label='Time ', valmin=0, valmax=115728, 
             valinit=0, valfmt='%i Frames', facecolor='#cc7000')
Play=Button(ax=ax_play, label='Play', image=None, color='0.85', hovercolor='0.95')
# Plot default data

im = ax.imshow(raw_data['data'][0])


# Update values
def update(val):
    T = s_T.val
    im.set_data(raw_data['data'][T])
    fig.canvas.draw_idle()
    
def mov(val):
    
        im.set_data(raw_data['data'][0])
        im.set_data(raw_data['data'][i])
        return (im,)
    
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=30, interval=1000./30, blit=True)


Play.on_clicked(mov)    
s_T.on_changed(update)