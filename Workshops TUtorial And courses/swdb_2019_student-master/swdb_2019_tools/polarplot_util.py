"""
Useful polar plotting functions
""" 
import matplotlib.pyplot as plt
import numpy as np

def polar_plot(ax, theta, r, color='black', linewidth = 2):  
  """function for making a polar plot to summarize a cell's directional preferences
    
    Parameters
    ----------
    ax : pyplot axis
      the axis onto which the polar plot will be drawn
    theta : list<Number>
      list of the directional preferences, in degrees.
    r : list<Number>
      list of preference magnitude. Range: [0,1].
      the number of elements should match 'theta'
    color : string
      the color of the polar plot line. E.g. 'black', 'grey'
    linewidth : Number
      the width of the polar plot line.
    
    Example
    -------
    theta = np.radians([0,90, 180, 270])
    r = [1.0,0.3, 0.5, 1.0]
    theta2 = np.radians([0,90, 180, 270])
    r2 = [0.5,0.5, 0.3, 0.5]
    theta3 = np.radians([0,90, 180, 270])
    r3 = [0.4,0.4, 0.4, 0.4]
        
    fig, axes = plt.subplots(nrows=1, ncols=1, subplot_kw=dict(polar=True))
    polar_plot(axes, theta2, r2, color='grey', linewidth=2)
    polar_plot(axes, theta3, r3, color='grey', linewidth=2)
    polar_plot(axes, theta, r, color='black', linewidth=4)
  """
  # Need to re-add the first point so we have an enclosed polygon.
  theta = np.append(theta, theta[0])
  r = np.append(r, r[0])
  ax.plot(theta, r, color=color, ls='-', linewidth=linewidth)

  # Because default is 0 pointing to east, and we want 0 at north.
  ax.set_theta_zero_location('N')
  # Default is ccw, we want cw
  ax.set_theta_direction(-1)

  ax.set_xticks(np.radians([0,45,90,135,180,225,270,315]))
  # How often to show the 0 ... 1.0 gradation in r values.
  # Below is a hack. For some reason without setting the y ticks,
  # the plots will get disfigured when you add more plots in the same axis.
  ax.set_yticks([0,1.0])
  ax.get_yaxis().set_visible(False)
  # Turn off the gray scaling circles and r grids
  ax.grid(False)
  
  # rmax has to be set after plotting. See https://stackoverflow.com/questions/54653423/matplotlib-set-rmax-and-set-rticks-not-working
  ax.set_rmax(1.0)