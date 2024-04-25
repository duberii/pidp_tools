import plotly.graph_objects as go 
import plotly.express as pxp
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from ipywidgets import interactive, HBox, VBox, Layout
import ipywidgets as widgets
import subprocess
from pidp_tools.analysis import get_charge
import math
from PIL import Image

def plot_vector(x_val,y_val,z_val):
  """
  Plots the vector with the supplied components.

  Parameters
  ----------
  x_val \: float
      The x component of the vector to be plotted.
  y_val \: float
      The y component of the vector to be plotted.
  z_val \: float
      The z component of the vector to be plotted.

  """
  x = [0, x_val]
  y = [0, y_val]
  z = [0, z_val]
  data = [go.Scatter3d(x=x, y=y, z=z,mode="lines",type="scatter3d",hoverinfo="none",line=dict(color="blue",width=3),name="Vector"),go.Cone(x=[x[1]],y=[y[1]],z=[z[1]],u=[0.3 * (x[1] - x[0])],v=[0.3 * (y[1] - y[0])],w=[0.3 * (z[1] - z[0])],anchor="tip",hoverinfo="none",colorscale=[[0, "blue"], [1, "blue"]],showscale=False),go.Scatter3d(x=[0],y=[0],z=[0],mode="markers",marker=dict(size=5,color="red"),name="Origin")]
  layout = go.Layout(scene=dict(aspectmode="cube",xaxis=dict(range=[-10,10]),yaxis=dict(range=[-10,10]),zaxis=dict(range=[-10,10])))
  fig = go.Figure(data=data, layout=layout)
  fig.show()

def x_coord(r_curvature, theta_curvature,center_x, t,q):
  """
  Calculates the x coordinate of a track based on its radius of curvature, starting point, center of curvature, and charge after t radians.

  Parameters
  ----------
  r_curvature \: float
      The radius of curvature of the track.
  theta_curvature \: float
      The phase offset of the track's parametrization. This is necessary to ensure the track starts at the origin.
  center_x \: float
      The x coordinate of the center of curvature of the track of the particle.
  t \: :external:func:`numpy.array`
      An array of angles along the track (in radians) to calculate the x coordinate of.
  q \: int
      The charge of the particle that left the track.
        
  Returns
  -------
  :external:func:`numpy.array`
      The x coordinates of the points along the track with angles provided in the t argument.
  """
  return r_curvature*np.cos(theta_curvature-q*t) + center_x

def y_coord(r_curvature,theta_curvature,center_y, t,q):
  """
  Calculates the y coordinate of a track based on its radius of curvature, starting point, center of curvature, and charge.

  Parameters
  ----------
  r_curvature \: float
      The radius of curvature of the track.
  theta_curvature \: float
      The phase offset of the track's parametrization. This is necessary to ensure the track starts at the origin.
  center_y \: float
      The y coordinate of the center of curvature of the track of the particle.
  t \: :external:func:`numpy.array`
      An array of angles along the track (in radians) to calculate the y coordinate of.
  q \: int
      The charge of the particle that left the track.
        
  Returns
  -------
  :external:func:`numpy.array`
      The y coordinates of the points along the track with angles provided in the t argument.
  """
  return r_curvature*np.sin(theta_curvature-q*t)+ center_y

def degrees_to_radians(degrees):
  """
  Converts degrees to radians.

  Parameters
  ----------
  degrees \: float or :external:func:`numpy.array`
      Angle (or array of angles) in degrees.
        
  Returns
  -------
  float or :external:func:`numpy.array`
      Angle (or array of angles) in radians.

  """
  return degrees*np.pi/180

class wirePositions():
  """
  An object that stores a variety of properties of the straws of the GlueX Central Drift Chamber (CDC).

  Attributes
  -----------
  wire_positions_df \: :external:class:`pandas.DataFrame`
      A pandas dataframe with at least the following three columns:

      - "pos" \: contains a list with two elements, each of which is a position vector. The first position vector is the position of the straw at the upstream end of the detector (z=s=0). The second position vector is the position of the straw at the downstream end of the detector (z=175, s=1). Positions in between can be expressed as a convex combination of the two vectors based on the normalized z coordinate.
      - "ring" \: The ring the straw in question belongs to
      - "straw" \: The straw in question
  positionsMatrix \: list
      A 3D jagged array of wire positions. The first index of the array represents the ring, and the second index represents the straw number. The third index represents the end of the detector (index 0 is upstream, index 1 is downstream). The fourth index represents the coordinate of the position of the straw in question (index 0 is x, index 1 is y, and index 2 is z).
  """
  def __init__(self):
    subprocess.run(["wget", "-q", "-O", "CentralDC_HDDS.xml", "https://github.com/JeffersonLab/hdds/raw/master/CentralDC_HDDS.xml"])
    tree =ET.parse("CentralDC_HDDS.xml")
    root = tree.getroot()
    n_wires_per_ring = [0]
    positionsMatrix = []
    ar = []
    for tag in [j for j in root.findall("composition") if j.get("name")=="CDClayers"][0]:
      if tag.tag == "mposPhi":
        positionsMatrix.append([])
        n_wires_per_ring.append(int(tag.get('ncopy')))
        R_Z = float(tag.get("R_Z").split()[0])
        Phi0 = degrees_to_radians(float(tag.get("Phi0")))
        if tag.get("volume") == 'CDCstrawShort':
          dPhi = degrees_to_radians(float(tag.get("dPhi")))
          for i in range(int(tag.get("ncopy"))):
            ar.append({'pos':[[R_Z*np.cos(Phi0+i*dPhi),R_Z*np.sin(Phi0+i*dPhi),0],[R_Z*np.cos(Phi0+i*dPhi),R_Z*np.sin(Phi0+i*dPhi),150]],'ring':int(tag.find("ring").get("value")),'straw':i+1})
            positionsMatrix[-1].append(ar[-1]['pos'])
        if tag.get("volume") == 'CDCstrawLong':
          n_wires_per_ring.append(int(tag.get('ncopy')))
          dPhi = 2*np.pi/int(tag.get("ncopy"))
          rot = degrees_to_radians(float(tag.get('rot').split()[0]))
          for i in range(int(tag.get("ncopy"))):
            ar.append({'pos':[[R_Z*np.cos(Phi0+i*dPhi)+75*np.tan(rot)*np.sin(Phi0+i*dPhi),R_Z*np.sin(Phi0+i*dPhi)-75*np.tan(rot)*np.cos(Phi0+i*dPhi),0],[R_Z*np.cos(Phi0+i*dPhi)-75*np.tan(rot)*np.sin(Phi0+i*dPhi),R_Z*np.sin(Phi0+i*dPhi)+75*np.tan(rot)*np.cos(Phi0+i*dPhi),150]],'ring':int(tag.find("ring").get("value")),'straw':i+1})
            positionsMatrix[-1].append(ar[-1]['pos'])
    self.wire_positions_df = pd.DataFrame.from_records(ar)
    self.positionsMatrix = positionsMatrix
  def position(self, rings, wires,s, dataframe=False):
    """
    Calculates the wire position.

    Parameters
    ------------
    rings \: list
        A list of rings that represent the rings in which the wires provided are located.
    wires \: list
        A list of wires whose positions will be calculated.
    s \: float or list
        A number (or list of numbers) that represents a normalized z coordinate (z/175) of the wire. This is used because several wires are slanted, so their x and y coordinates depend on the z coordinate. If s is a float or int, all wire positions will be calculated with this same normalized z coordinate.
    t \: :external:func:`numpy.array`
      An array of angles along the track (in radians) to calculate the x coordinate of
    dataframe \: bool, default False
      Return a dataframe containing the x and y coordinates of the wires. If True, the coordinates will be returned in a dataframe with four columns\: "ring", "wire", "x", and "y".
        
    Returns
    ---------
    :external:func:`numpy.array` or :external:class:`pandas.DataFrame`
        An array or dataframe containing the x and y coordinates of the wires at the given normalized z coordinate(s).
    """
    if dataframe:
      df = pd.DataFrame()
      df['ring'] = rings
      df['wire'] = wires
      if not isinstance(s,float) and not isinstance(s, int):
        df['x'] = [self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s[i]) + s[i]* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))]
        df['y'] = [self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s[i]) + s[i]* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))]
      else:
        df['x'] = [self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s) + s* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))]
        df['y'] = [self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s) + s* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))]
      return df
    else:
      if not isinstance(s,float) and not isinstance(s, int):
        new_x = np.array([self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s[i]) + s[i]* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))])
        new_y = np.array([self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s[i]) + s[i]* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))])
      else:
        new_x = np.array([self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s) + s* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))])
        new_y = np.array([self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s) + s* self.positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))])
      return new_x, new_y

class track_fitting():
  """
  Creates a track fitting applet with hits drawn based on the provided event.

  Parameters
  ----------
  title \: str, default "Hits in the GlueX CDC"
      The title of the plot.
  wire_positions \: :class:`wirePositions`, default None
      A wirePositions object to reference detector coordinates. By default, one will be created. However, if you wish to avoid creating multiple instances, you can provide one that is already created.
  event \: :external:class:`pandas.DataFrame`, default None
      A row of a dataframe with the following columns:

      - "px" \: The x component of the momentum of the generated particle.
      - "py" \: The y component of the momentum of the generated particle.
      - "pz" \: The z component of the momentum of the generated particle.
      - "vz" \: The z coordinate of the vertex ("creation point") of the generated particle.
      - "ring" \: A list of rings that contain hits.
      - "straw" \: A list of straws that contain hits.
  showTrack \: bool, default False
      Draws the track. If True, a track will be drawn in red.
  showlegend \: bool, default False
      Displays a legend. If True, a legend will be drawn to the right of the applet.
        
  Attributes
  ----------
  figure \: :external:class:`plotly.graph_objects.Figure`
      A plotly FigureWidget that contains the applet.
  px \: float
      The x component of the momentum of the generated particle.
  py \: float
      The y component of the momentum of the generated particle.
  pz \: float
      The z component of the momentum of the generated particle.
  current_px \: float
      The x component of the momentum, as determined by the current state of any available sliders.
  current_py \: float
      The y component of the momentum, as determined by the current state of any available sliders.
  current_pz \: float
      The z component of the momentum, as determined by the current state of any available sliders.
  current_z0 \: float
      The z coordinate of the vertex ("creation point"), as determined by the current state of any available sliders.
  current_charge \: int or float
      The charge, as determined by the current state of any available sliders.
  ring \: list
      A list of rings that contain hits.
  straw \: list
      A list of straws that contain hits.
  vz \: float
      The z coordinate of the vertex ("creation point") of the generated particle.
  wirePositions \: :class:`wirePositions`
      A wirePositions object that contains information about the GlueX CDC.

  """
  def __init__(self, title="Hits in the GlueX CDC", wire_positions=None, event=None, showTrack = False, showlegend=False):
    if wire_positions is None:
      self.wirePositions = wirePositions()
    else:
      self.wirePositions = wire_positions
    self.figure = go.FigureWidget()
    self.figure.update_layout(xaxis_range=[-60,60], yaxis_range=[-60,60],width=500,height=500,showlegend=showlegend,title=title,xaxis_title="X", yaxis_title="Y")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-10.5, y0=-10.5, x1=10.5, y1=10.5, line_color="black")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-55, y0=-55, x1=55, y1=55, line_color="black")
    if event is not None:
      self.charge = get_charge(event["particle"])
      self.px = event['px']
      self.py = event['py']
      self.pz = event['pz']
      self.vz = event['vz']
      self.ring = event['ring']
      self.straw = event['straw']
      df = self.wirePositions.position(self.ring,self.straw,75,dataframe=True)
      self.figure.add_scatter(x=df['x'].to_numpy(),y=df['y'].to_numpy(),mode='markers',text = ["Ring: " + str(row['ring']) + "\n Wire: " + str(row['wire']) for row in df.iloc],hoverinfo='text',marker={"size":3})
    else:
      self.px = 0.
      self.py = 0.
      self.pz = 0.
      self.vz = 0.
      self.ring = []
      self.straw = []
    self.current_px = 0.
    self.current_py = 0.
    self.current_pz = 0.
    self.current_z0 = 0.
    self.current_charge = 1
    if showTrack:
      self.figure.add_scatter(mode='lines',name='track',marker={"color":'red'})
  def show(self):
    """
    Shows the track fitting applet. This is a wrapper method for `plotly.graph_objects.Figure.show <https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html#plotly.graph_objects.Figure.show>`_.
    
    """
    self.figure.show()
  def update_layout(self, *args, **kwargs):
    """
    Updates the layout of the applet. This is a wrapper method for :external:meth:`plotly.graph_objects.Figure.update_layout`. Refer to its documentation for a list of all possible arguments.
    
    """
    self.figure.update_layout(*args, **kwargs)
  def update(self, *args, **kwargs):
    """
    Updates the applet. This is a wrapper method for :external:meth:`plotly.graph_objects.Figure.update`. Refer to its documentation for a list of all possible arguments.
    
    """
    self.figure.update(*args, **kwargs)
  def display(self,**kwargs):
    """
    Displays the sliders corresponding to track parameters. This is done by calling the :meth:`track_fitting.update_figure` method.

    Parameters:
    -----------
    \*\*kwargs:
        The following keyword arguments correspond to the arguments of the :meth:`track_fitting.update_figure` method. Each of these keyword arguments should be used for every call. Each of these keyword arguments can be assigned one of the following:

        - A tuple of the form `(min_value, max_value, step)` for numerical (float/int) arguments, which will generate a slider for the argument in question.
        - A boolean for boolean arguments, which will generate a checkbox, which allows the user to toggle the argument.
        - A list of strings or tuples, which will produce a dropdown menu of the provided options. If strings, the values in the dropdown will be passed directly to the :meth:`track_fitting.update_figure` method. If a list of tuples is passed, they should be of the form (label, value), where the dropdown options are the "labels", and the corresponding value to be passed to :meth:`track_fitting.update_figure` is the "value".
        - An instance of `ipywidgets.widgets.interaction.fixed <https://ipywidgets.readthedocs.io/en/stable/reference/ipywidgets.html#ipywidgets.widgets.interaction.fixed>`_, which represents that the provided value cannot be changed by the user. No widget will be generated for this argument.
          
        The below list of keyword arguments specifies the expected type of the arguments.
        
            px \: int or float
                The x component of the momentum of the track.
            py \: int or float
                The y component of the momentum of the track.
            pz \: int or float
                The z component of the momentum of the track.
            z0 \: int or float
                The z coordinate of the vertex that produced the track. In other words, this is the z component of the "starting point" of the track.
            charge \: int or float
                The charge of the track. This can be a float, but it should usually be -1 or 1.
            t_max \: int or float
                The length of the track to draw, measured in radians.
            show_track \: bool
                Draws the track. If True, draws the track with the supplied parameters.
            show_hits \: bool
                Draws the hits. If true, draws the hits detailed in the event passed during initialization.
    
    Examples:
    -----------
    >>> fig = track_fitting()
    >>> fig.update()
    >>> fig.display(px=(-1,1,0.0001),py=(-1,1,0.0001),t_max=(0,2*np.pi,0.0001),charge = (-1,1,2),show_track=fixed(True),show_hits=fixed(False),pz=fixed(0),z0=fixed(75))

    The result of this code is a track fitting applet with 4 sliders\: px, py, t_max, and charge. The remaining arguments are fixed, and will not correspond to a displayed widget.
    
    """
    widget=interactive(self.figure.update_figure, **kwargs)
    controls = HBox(widget.children[:-1], layout = Layout(flex_flow='row wrap'))
    output = widget.children[-1]
    display(VBox([controls, output]))
  def update_figure(self, z0=0,px=0,py=0,charge=1,pz=0,t_max = np.pi, show_track = True, show_hits = True):
    """
    Updates and draws the figure, recalculating the positions of the hits and the track trajectory, based on the provided track parameters.

    Parameters:
    -----------
    px \: float, default 0
        The x component of the momentum of the track.
    py \: float, default 0
        The y component of the momentum of the track.
    pz \: float, default 0
        The z component of the momentum of the track.
    z0 \: float, default 0
        The z coordinate of the vertex that produced the track. In other words, this is the z component of the "starting point" of the track.
    charge \: float or int, default 1
        The charge of the track. This can be a float, but it should usually be -1 or 1.
    t_max \: float, default 3.14
        The length of the track to draw, measured in radians.
    show_track \: bool, default True
        Draws the track. If True, draws the track with the supplied parameters.
    show_hits \: bool, default True
        Draws the hits. If true, draws the hits detailed in the event passed during initialization.

    """
    self.current_px = px
    self.current_py = py
    self.current_pz = pz
    self.current_z0 = z0
    self.current_charge = charge

    p = (self.current_px**2 + self.current_py**2 + self.current_pz**2)**0.5
    try:
      z_angle = np.arcsin(pz/p)
    except:
      z_angle = 0
    center_x = self.current_charge*330*self.current_py/1.7
    center_y = -1*self.current_charge*330*self.current_px/1.7
    r_curvature = (center_x**2 + center_y**2)**0.5
    theta_curvature = np.arctan2(-1*center_y,-1*center_x)
    if show_hits:
      new_x,new_y = self.wirePositions.position(self.ring,self.straw,z0/175)
      theta_hits= np.arctan2(new_y-center_y,new_x-center_x)
      for i in range(len(theta_hits)):
        theta_hits[i]-= theta_curvature
        theta_hits[i] = theta_hits[i] * -1 * charge
        if theta_hits[i] <0:
          theta_hits[i] += 2*np.pi
      for i in range(len(self.ring)):
        new_x,new_y = self.wirePositions.position(self.ring,self.straw,[(z0 +r_curvature*theta_hits[i]*np.tan(z_angle))/175 for i in range(len(self.ring))])
    if show_track:
      t_curve = np.linspace(0,t_max,1000)
      x_track = x_coord(r_curvature,theta_curvature,center_x,t_curve,charge)
      y_track = y_coord(r_curvature,theta_curvature,center_y,t_curve,charge)
    with self.figure.batch_update():
      if show_track:
        self.figure.data[len(self.figure.data)-1]['x']=x_track
        self.figure.data[len(self.figure.data)-1]['y']=y_track
      self.update_layout(width=500,height=500,showlegend=False)
      if show_hits:
        self.figure.data[0]['x']=new_x
        self.figure.data[0]['y']=new_y
    self.show()

  def distance_to_curve(self,x,y):
    """
    Calculates the distance between a point and the "correct" track. Used to find RMSE by varying the z coordinate of the vertex.

    Parameters:
    -----------
    x \: float
        The x coordinate of the point.
    y \: float
        The y coordinate of the point.

    Returns:
    --------
    distance_to_nearest_curve_point \: float
        The distance between the supplied point and the "correct" track, as determined by the parameters provided by the event.
    
    """
    center_x = self.charge*330*self.py/1.7
    center_y = -330*self.charge*self.px/1.7
    r_curvature = (center_x**2 + center_y**2)**0.5
    distance_to_nearest_curve_point = abs(r_curvature-((x-center_x)**2+(y-center_y)**2)**0.5)
    return distance_to_nearest_curve_point
  
  def calculate_rmse(self,x_points, y_points):
    """
    Calculates the RMSE of the distance between the points and the "correct" track. This is used to optimize the z coordinate of the vertex.

    Parameters:
    -----------
    x_points \: list
        A list of x coordinates of hits in the CDC.
    y_points \: list
        A list of y coordinates of hits in the CDC.

    Returns:
    --------
    rmse \: float
        The root mean square error (RMSE) of the distances of the hits to the "correct" track.

    """
    mse = sum([self.distance_to_curve(x_points[j],y_points[j])**2 for j in range(len(x_points))])/len(x_points)
    rmse = mse**0.5
    return rmse
  
  def calc_z0(self, quiet=False):
    """
    Finds the most likely z coordinate of the vertex by minimizing the RMSE of the hit positions.

    Parameters:
    -----------
    quiet \: bool, default False
        Suppresses printing of z coordinate. If True, the correct z0 is simply returned and not printed.

    Returns:
    --------
    rmse \: float
        The root mean square error (RMSE) of the distances of the hits to the "correct" track. Only returned if `quiet` is True.
    
    """
    rmses = {}
    for z in np.linspace(0,150,151):
      new_x, new_y = self.wirePositions.position(self.ring, self.straw,z/175)
      rmses[z] = self.calculate_rmse(new_x,new_y)
    if quiet:
      return min(rmses, key=rmses.get)
    else:
      print("The true z value is about " + str(min(rmses, key=rmses.get)))

  def show_answer(self, charge=True,px=True,py=True,pz=False,z0=False):
    """
    Prints the track parameters, as determined by the GlueX reconstruction algorithm. The z coordinate of the vertex is the one exception, as this is calculated by minimizing RMSE.

    Parameters:
    -----------
    charge \: bool, default True
        Prints the charge of the particle. If True, the charge of the particle is printed
    px \: bool, default True
        Prints the x component of the momentum of the track. If True, the x component of momentum is printed.
    py \: bool, default True
        Prints the y component of the momentum of the track. If True, the y component of momentum is printed.
    pz \: bool, default False
        Prints the z component of the momentum of the track. If True, the y component of momentum is printed.
    z0 \: bool, default False
        Prints the z component of the vertex of the track. If True, the z component of the track vertex is printed.

    """
    if z0:
      self.calc_z0()
    if charge:
      print("The charge of the particle is " + str(self.charge))
    if px:
      print("The x component of momentum is " + str(round(self.px, 2)))
    if py :
      print("The y component of momentum is " + str(round(self.py, 2)))
    if pz:
      print("The z component of momentum is " + str(round(self.pz, 2)))
  def __repr__():
    return ""

class CDC_plot_2D():
  """
  Creates a 2D model of the straws of the CDC.

  Parameters:
  -----------
  rings \: list
      A list of rings to plot. All straws from the selected ring will be drawn. Each ring should be an integer between 1 and 28.
  title \: str, default "Wire Positions in the GlueX CDC"
      The title of the figure.
  wire_positions \: :meth:`wirePositions`, default None
      A :meth:`wirePositions` object to reference detector coordinates. By default, one will be created. However, if you wish to avoid creating multiple instances, you can provide one that is already created.
  showlegend \: bool, default True
      Displays a legend to the right of the model. If True, displays the legend.

  Attributes:
  -----------
  figure \: :external:class:`plotly.graph_objects.Figure`
      A plotly figure that displays the requested rings of the GlueX CDC.

  """
  def __init__(self,rings,title="Wire Positions in the GlueX CDC", wire_positions=None,showlegend=True):
    if wire_positions is None:
      self.wirePositions = wirePositions()
    else:
      self.wirePositions = wire_positions
    self.figure = go.FigureWidget()
    self.figure.update_layout(xaxis_range=[-60,60], yaxis_range=[-60,60],width=500,height=500,showlegend=showlegend,title=title,xaxis_title="X", yaxis_title="Y")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-10.5, y0=-10.5, x1=10.5, y1=10.5, line_color="black")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-55, y0=-55, x1=55, y1=55, line_color="black")
    for ring in rings:
      self.figure.add_scatter(mode='markers',name="Ring " + str(ring),marker={"size":3},hoverinfo='text',text = ['Ring: ' + str(ring) + "\n Wire: " + str(k) for k in range(len(self.wirePositions.positionsMatrix[int(ring)-1]))])
  
  def show(self):
    """
    Shows the CDC plot track fitting applet. This is a wrapper method for `plotly.graph_objects.Figure.show <https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html#plotly.graph_objects.Figure.show>`_

    """
    self.figure.show()
  
  def update_layout(self, *args, **kwargs):
    """
    Updates the layout of the applet. This is a wrapper method for :external:meth:`plotly.graph_objects.Figure.update_layout`. Refer to its documentation for a list of all possible arguments.
    
    """
    self.figure.update_layout(*args, **kwargs)
  
  def update(self, *args, **kwargs):
    """
    Updates the applet. This is a wrapper method for :external:meth:`plotly.graph_objects.Figure.update`. Refer to its documentation for a list of all possible arguments.
    
    """
    self.figure.update(*args, **kwargs)
  def display(self,**kwargs):
    """
    Displays the 2D CDC plot. This is done by calling the :meth:`CDC_plot_2D.update_figure` method.

    Parameters:
    -----------
    \*\*kwargs:
        The following keyword arguments correspond to the arguments of the :meth:`track_fitting.update_figure` method. Each of these keyword arguments should be used for every call. Each of these keyword arguments can be assigned one of the following:

        - A tuple of the form "(min_value, max_value, step)" for numerical (float/int) arguments, which will generate a slider for the argument in question.
        - A boolean for boolean arguments, which will generate a checkbox, which allows the user to toggle the argument.
        - A list of strings or tuples, which will produce a dropdown menu of the provided options. If strings, the values in the dropdown will be passed directly to the :meth:`CDC_plot_2D.update_figure` method. If a list of tuples is passed, they should be of the form (label, value), where the dropdown options are the "labels", and the corresponding value to be passed to :meth:`CDC_plot_2D.update_figure` is the "value".
        - An instance of `ipywidgets.widgets.interaction.fixed <https://ipywidgets.readthedocs.io/en/stable/reference/ipywidgets.html#ipywidgets.widgets.interaction.fixed>`_, which represents that the provided value cannot be changed by the user. No widget will be generated for this argument.
  
        The below list of keyword arguments specifies the expected type of the arguments.
        
        z \: float or int
            The z coordinate of the slice of the CDC you want to look at. Should be between 0 and 175 (this quantity is assumed to be in cm).
        rings \: list
            A list of rings to be displayed. Each of the elements of this list should be an int, and should be between 1 and 28.
    
    Examples:
    -----------
    >>> fig = CDC_plot_2D()
    >>> fig.update()
    >>> fig.display(z=(0,175,0.1), rings=fixed([1,2,3]))

    The result of this code is a plot with a slider that allows the user to look at different z-slices of the CDC. The rings shown are fixed.
    
    """
    widget=interactive(self.update_figure, **kwargs)
    controls = HBox(widget.children[:-1], layout = Layout(flex_flow='row wrap'))
    output = widget.children[-1]
    display(VBox([controls, output]))

  def update_figure(self,z=0, rings=[]):
    """
    Updates and draws the figure, recalculating the positions of the straws based on the provided z coordinate.

    Parameters:
    -----------
    z \: numerical
        The z coordinate of the slice of the CDC you want to look at. Should be between 0 and 175 (this quantity is assumed to be in cm).
    rings \: list
        A list of rings to be displayed. Each of the elements of this list should be an int, and should be between 1 and 28.
    
    """
    s= z/175
    with self.figure.batch_update():
      for i in range(len(rings)):
        self.figure.data[i]['x']=[(1-s)*wire[0][0]+ s*wire[1][0] for wire in self.wirePositions.positionsMatrix[int(rings[i])-1]]
        self.figure.data[i]['y']=[(1-s)*wire[0][1]+ s*wire[1][1] for wire in self.wirePositions.positionsMatrix[int(rings[i])-1]]
    self.show()

  def __repr__(self):
    return ""

class CDC_plot_3D():
  """
  Creates and displays a 3D model of the straws of the CDC.

  Parameters:
  -----------
  rings \: list
      A list of rings to plot. All straws from the selected ring will be drawn. Each ring should be an integer between 1 and 28. Can supply at most 3 rings.
  wire_positions \: :class:`wirePositions`, default None
      A wirePositions object to reference detector coordinates. By default, one will be created. However, if you wish to avoid creating multiple instances, you can provide one that is already created.

  Attributes:
  -----------
  figure \: :external:class:`plotly.graph_objects.Figure`
      A plotly figure that displays the requested rings of the GlueX CDC.

  """
  def __init__(self, rings, wire_positions = None):
    if wire_positions is None:
      self.wirePositions = wirePositions()
    else:
      self.wirePositions = wire_positions
    if len(rings) > 3:
      raise ValueError("Too many rings. Provide at most 3 rings.")
    xs = []
    ys = []
    zs = []
    ring_wires = []
    rings_to_plot = []
    for i in self.wirePositions.wire_positions_df.iloc:
      if i['ring'] in rings:
        xs.append(i['pos'][0][0])
        xs.append(i['pos'][1][0])
        ys.append(i['pos'][0][1])
        ys.append(i['pos'][1][1])
        zs.append(i['pos'][0][2])
        zs.append(i['pos'][1][2])
        ring_wires.append("Ring: "+ str(i['ring']) + " Wire: " + str(i['straw']))
        ring_wires.append("Ring: "+ str(i['ring']) + " Wire: " + str(i['straw']))
        rings_to_plot.append(i['ring'])
        rings_to_plot.append(i['ring'])
    data_to_plot = pd.DataFrame()
    data_to_plot['x'] = xs
    data_to_plot['y'] = ys
    data_to_plot['z'] = zs
    data_to_plot['Ring and Wire']=ring_wires
    data_to_plot['ring']=rings_to_plot
    self.figure = pxp.line_3d(data_frame=data_to_plot,x='x',y='z',z='y',line_group='Ring and Wire',color='ring')
    self.figure.update_traces(line=dict(width=5))
    self.figure.update_layout(title="Rings of the GlueX Central Drift Chamber",legend_title="Ring Number",scene={'aspectmode':'cube','xaxis':{'range':[-60,60],"title":"X"},'yaxis':{'range':[0,175],"title":"Z"},'zaxis':{'range':[-60,60],'title':'Y'}})
    self.figure.show()

  def __repr__(self):
    return ""

class interactive_image():
  """
  Creates an interactive image that can be used to find exponential cuts in 2D histograms.

  Parameters
  ----------
  image_path \: str, default ""
      The path to the saved ROOT canvas.
  x_max \: float, default 1.
      The maximum of the x axis, which is typically used for momentum. This should be in the same units as were used when creating the ROOT histogram.
  y_max \: float, default 2.*10**-5
      The maximum of the y_axis. This should be in the same units as were used when creating the ROOT histogram.
        
  Attributes
  ----------
  figure \: :external:class:`plotly.graph_objects.Figure`
      A plotly FigureWidget that contains the interactive image.
  x_max \: float, default 1.
      The maximum of the x axis, which is typically used for momentum. This should be in the same units as were used when creating the ROOT histogram.
  y_max \: float, default 2.*10**-5
      The maximum of the y_axis. This should be in the same units as were used when creating the ROOT histogram.

  """
  def __init__(self,image_path="",x_max = 1., y_max=2.*10**-5):
    self.figure = go.FigureWidget()
    self.x_max = x_max
    self.y_max = y_max
    self.figure.add_scatter(x=np.linspace(0,self.x_max,1000), y=np.zeros(1000),mode='lines',name='Cut',marker={"color":'red'})
    print("This equation is rescaled by 10**" + str(math.floor(math.log(self.y_max,10))-1))
    if len(image_path)>0:
      self.figure.add_layout_image({"source":Image.open(image_path),"xref":"x", "yref":"y","x": -0.125*self.x_max,"y":1.125*self.y_max, "sizex": 1.25*self.x_max,"sizey":1.25*self.y_max,"sizing":"stretch","opacity":1,"layer":"below"})
    self.figure.update_xaxes(range=[0, self.x_max],showgrid=False)
    self.figure.update_yaxes(range=[0, self.y_max],showgrid=False)
    self.figure.update_layout(template="plotly_white",showlegend=False, autosize=False,width=530,height=457,margin=dict(l=30,r=0,b=50,t=50,pad=0), title="Ionization Energy Loss vs. Momentum", xaxis_title="Momentum (GeV/c)", yaxis_title="dE/dx")
  
  def show(self):
    """
    Shows the interactive image. This is a wrapper method for `plotly.graph_objects.Figure.show <https://plotly.com/python-api-reference/generated/plotly.graph_objects.Figure.html#plotly.graph_objects.Figure.show>`_

    """
    self.figure.show()
  def update_layout(self, *args, **kwargs):
    """
    Updates the layout of the applet. This is a wrapper method for :external:meth:`plotly.graph_objects.Figure.update_layout`. Refer to its documentation for a list of all possible arguments.
    
    """
    self.figure.update_layout(*args, **kwargs)
  def update(self, *args, **kwargs):
    """
    Updates the layout of the applet. This is a wrapper method for :external:meth:`plotly.graph_objects.Figure.update_layout`. Refer to its documentation for a list of all possible arguments.
    
    """
    self.figure.update(*args, **kwargs)
  
  def display(self, **kwargs):
    """
    Displays the sliders for the interactive image. This is done by calling the :meth:`interactive_image.update_figure` method.

    Parameters:
    -----------
    \*\*kwargs:
        The following keyword arguments correspond to the arguments of the :meth:`interactive_image.update_figure` method. Each of these keyword arguments should be used for every call. Each of these keyword arguments can be assigned one of the following:
              
        - A tuple of the form `(min_value, max_value, step)` for numerical (float/int) arguments, which will generate a slider for the argument in question.
        - A boolean for boolean arguments, which will generate a checkbox, which allows the user to toggle the argument.
        - A list of strings or tuples, which will produce a dropdown menu of the provided options. If strings, the values in the dropdown will be passed directly to the :meth:`interactive_image.update_figure` method. If a list of tuples is passed, they should be of the form (label, value), where the dropdown options are the "labels", and the corresponding value to be passed to :meth:`interactive_image.update_figure` is the "value".
        - An instance of `ipywidgets.widgets.interaction.fixed <https://ipywidgets.readthedocs.io/en/stable/reference/ipywidgets.html#ipywidgets.widgets.interaction.fixed>`_, which represents that the provided value cannot be changed by the user. No widget will be generated for this argument.
            
        The below list of keyword arguments specifies the expected type of the arguments. Note that these correspond to the constants in the equation y= exp(a*p+b)+c. Also note that the y axis is automatically rescaled, and the rescaling factor will be printed out when initializing the image.
        
        a \: numerical
            The exponential decay constant of the exponential function (e**(-a*p)).
        b \: numerical
            A vertical rescaling factor for the exponential function (by a factor of e**b).
        c \: numerical
            The vertical offset of the exponential function.

    Examples:
    -----------
    >>> fig = interactive_image()
    >>> fig.update()
    >>> fig.display(a=(-10,-1, 0.1), b=(0,10,0.1), c=fixed(1))

    The result of this code is a plot with a slider that allows the user to modify the values of a and b. The value of c will be fixed, so no slider will be displayed.
    
    """
    widget=interactive(self.update_figure, **kwargs)
    controls = HBox(widget.children[:-1], layout = Layout(flex_flow='row wrap'))
    output = widget.children[-1]
    display(VBox([controls, output]))

  def update_figure(self,a=1, b=1, c=1):
    """
    Updates and draws the interactive image, recalculating exponential equation y= exp(a*p+b)+c.

    Parameters:
    -----------
    a \: numerical
        The exponential decay constant of the exponential function (e**(-a*p)).
    b \: numerical
        A vertical rescaling factor for the exponential function (by a factor of e**b).
    c \: numerical
        The vertical offset of the exponential function.
    
    """
    with self.figure.batch_update():
      self.figure.data[0]['y']=10**(math.floor(math.log(self.y_max,10))-1)*(np.exp(a*np.linspace(0,self.x_max,1000)+b)+c)
    self.show()


def BCal_plot_2D():
  """
  Displays a 2D plot of the modules of the Barrel Calorimeter (BCal).

  """
  module_points = {i:[] for i in range(48)}
  inner_radius = 65
  outer_radius = 87.46
  inner_radius = 65
  inner_module_width = 8.52
  outer_module_width = 11.46
  fig = go.FigureWidget()
  fig.update_layout(xaxis_range=[-100,100], yaxis_range=[-100,100],width=520,height=500,xaxis_title="X", yaxis_title="Y")
  for i in range(48):
    angle = i * -1*np.pi/24 + np.pi
    unit_vector = np.array([np.cos(angle+np.pi/2),np.sin(angle+np.pi/2)])
    inner_radius_vector = inner_radius * np.array([np.cos(angle),np.sin(angle)])
    outer_radius_vector = outer_radius * np.array([np.cos(angle),np.sin(angle)])
    module_points[i].append(inner_radius_vector + inner_module_width/2 * unit_vector)
    module_points[i].append(inner_radius_vector - inner_module_width/2 * unit_vector)
    module_points[i].append(outer_radius_vector - outer_module_width/2 * unit_vector)
    module_points[i].append(outer_radius_vector + outer_module_width/2 * unit_vector)
    module_points[i].append(inner_radius_vector + inner_module_width/2 * unit_vector)
    fig.add_scatter(x=[module_points[i][j][0] for j in range(5)],y=[module_points[i][j][1] for j in range(5)],fill="toself",name="Module " + str(i+1),mode="lines")
  fig.show()

def BCal_plot_3D():
  """
  Displays a 3D model of the modules of the Barrel Calorimeter (BCal).

  """
  module_points = {i:[] for i in range(48)}
  inner_radius = 65
  outer_radius = inner_radius + 22.46
  inner_module_width = 8.52
  outer_module_width = 11.46
  fig = go.FigureWidget()
  fig.update_layout(scene = {"xaxis_title":'Z', "yaxis_title":'X', 'zaxis_title':'Y'})
  for module in range(48):
    angle = module * np.pi/24
    unit_vector = np.array([np.cos(angle+np.pi/2),np.sin(angle+np.pi/2),0])
    inner_radius_vector = inner_radius * np.array([np.cos(angle),np.sin(angle),0])
    outer_radius_vector = outer_radius * np.array([np.cos(angle),np.sin(angle),0])
    module_points[module].append(inner_radius_vector + inner_module_width/2 * unit_vector)
    module_points[module].append(inner_radius_vector - inner_module_width/2 * unit_vector)
    module_points[module].append(outer_radius_vector - outer_module_width/2 * unit_vector)
    module_points[module].append(outer_radius_vector + outer_module_width/2 * unit_vector)
    module_points[module].append(inner_radius_vector + inner_module_width/2 * unit_vector + np.array([0,0,390]))
    module_points[module].append(inner_radius_vector - inner_module_width/2 * unit_vector + np.array([0,0,390]))
    module_points[module].append(outer_radius_vector - outer_module_width/2 * unit_vector + np.array([0,0,390]))
    module_points[module].append(outer_radius_vector + outer_module_width/2 * unit_vector + np.array([0,0,390]))
    i = [0,0,0,3,0,1,1,2,2,3,4,4]
    j = [1,2,4,4,1,4,2,5,3,6,5,6]
    k = [2,3,3,7,4,5,5,6,6,7,6,7]
    fig.add_trace(go.Mesh3d(y=[pt[0] for pt in module_points[module]],z=[pt[1] for pt in module_points[module]],x=[pt[2] for pt in module_points[module]],name="Module " + str(module+1),flatshading=True,i=i,j=j,k=k,hovertemplate="Module " + str(module+1) + "<extra></extra>"))
  fig.show()

def BCal_module():
  """
  Displays a 2D model of a single Barrel Calorimeter (BCal) module, showing sector and layer segmentations.

  """
  def calculate_x(layer,sector):
    return 8.51*sector/4 + 1.475*(sector-2)/2*((layer+1)**2-layer-1)/20
  def calculate_y(layer,sector):
    return 22.46*((layer+1)**2-layer-1)/20
  fig = go.FigureWidget()
  fig.update_layout(xaxis_range=[-3,10], yaxis_range=[-1,25],width=600,height=600,xaxis_title="X", yaxis_title="Y",showlegend=False)
  fig.update_yaxes(scaleanchor="x",scaleratio=1)
  for layer in range(1,5):
    for sector in range(1,5):
      fig.add_scatter(x=[calculate_x(layer,sector),calculate_x(layer,sector-1),calculate_x(layer-1,sector-1),calculate_x(layer-1,sector)],y=[calculate_y(layer,sector),calculate_y(layer,sector-1),calculate_y(layer-1,sector-1),calculate_y(layer-1,sector)],fill="toself", mode="lines",name="Layer: " + str(layer) + ", Sector: " + str(sector),hoverinfo="text")
      fig.update()
  fig.show()


def view_shower_2D(df, eventNum):
  """
  Creates a displays a 2D visualization of the BCal shower in the provided row of the provided dataframe.

  Parameters:
  -----------
  df \: :external:class:`pandas.DataFrame`
      The dataframe containing events in the GlueX BCal. This dataframe should have 7 columns: "module","sector", "layer", "pulse_integral_up", "pulse_integral_down", "t_up", and "t_down". Each of these columns should contain a list, and all of these lists should be the same length across columns. Each element of each list corresponds to a single hit in the BCal.
  eventNum \: int
      The row of the dataframe for which you want to plot the shower.

  """
  def closest_hit(prev_layer,x,y):
    prev_distance = 10000000
    solution = (0,0)
    for hit in prev_layer.iloc:
      distance = ((x-hit['x'])**2 + (y-hit['y'])**2)**0.5
      if distance < prev_distance:
        solution = (hit['x'],hit['y'])
        prev_distance = distance
    if prev_distance > ((x)**2 + (y)**2)**0.5:
      return (x,y)
    else:
      return solution
  event = df.iloc[eventNum]
  inner_radius = 65
  outer_radius = 87.46
  inner_module_width = 8.52
  fig = go.FigureWidget()
  fig.update_layout(xaxis_range=[-100,100], yaxis_range=[-100,100],width=520,height=500,xaxis_title="X", yaxis_title="Y")
  points = {'x':[],'y':[],'z':[]}
  energies = []
  for i in range(len(event['sector'])):
    angle = event['module'][i] * -1*np.pi/24 + np.pi
    unit_radius_vector = np.array([np.cos(angle),np.sin(angle)])
    unit_vector = np.array([np.cos(angle+np.pi/2),np.sin(angle+np.pi/2)])
    inner_radius_vector = inner_radius * np.array([np.cos(angle),np.sin(angle)])
    origin = inner_radius_vector - inner_module_width/2 * unit_vector
    position = origin + unit_vector * (8.51*(event['sector'][i]-0.5)/4 + 1.475*(event['sector'][i]-0.5)/2*((event['layer'][i]+0.5)**2-event['layer'][i]-0.5)/20)
    position += 22.46*((event['layer'][i]+0.5)**2-event['layer'][i]-0.5)/20 *unit_radius_vector
    points['x'].append(position[0])
    points['y'].append(position[1])
    z = 1/30*16.75*(event['t_up'][i]-event['t_down'][i])
    points['z'].append(z)
    energies.append(2/130000*(event['pulse_integral_up'][i]*np.exp(z/525)+event['pulse_integral_down'][i]*np.exp((390-z)/525)))
  each_hit_df = pd.DataFrame()
  each_hit_df['module'] = event['module']
  each_hit_df['sector'] = event['sector']
  each_hit_df['layer'] = event['layer']
  each_hit_df['x'] = points['x']
  each_hit_df['y'] = points['y']
  each_hit_df['z'] = points['z']
  each_hit_df['E'] = energies
  each_hit_df['opacity'] = 0.75*(np.array(energies)-min(energies))/(max(energies)-min(energies))+0.25
  opacities = np.array(energies)/max(energies)
  previous_layer_hits = pd.DataFrame()
  for i in range(1,5):
    layer_hits = each_hit_df.loc[each_hit_df['layer']==i]
    if i != 1 and len(previous_layer_hits.index)>0 and len(layer_hits.index)>0:
      for hit in layer_hits.iloc:
        angle = hit['module'] * -1*np.pi/24 + np.pi
        unit_radius_vector = np.array([np.cos(angle),np.sin(angle)])
        unit_vector = np.array([np.cos(angle+np.pi/2),np.sin(angle+np.pi/2)])
        inner_radius_vector = inner_radius * np.array([np.cos(angle),np.sin(angle)])
        origin = inner_radius_vector - inner_module_width/2 * unit_vector
        position = origin + unit_vector * (8.51*(hit['sector']-0.5)/4 + 1.475*(hit['sector']-0.5)/2*((hit['layer']+0.5)**2-hit['layer']-0.5)/20)
        position += 22.46*((hit['layer']+0.5)**2-hit['layer']-0.5)/20 *unit_radius_vector
        x2 = position[0]
        y2= position[1]
        x1,y1 = closest_hit(each_hit_df.loc[each_hit_df['layer']<i],x2,y2)
        fig.add_shape(type="line",x0=x1,y0=y1,x1=x2,y1=y2,opacity=hit['opacity'],line_color="#636EFA")
    previous_layer_hits = layer_hits
  fig.add_scatter(x=points['x'],y=points['y'],mode='markers',marker = {'opacity':opacities})
  fig.add_shape(type="circle", xref="x", yref="y", x0=-1*inner_radius, y0=-1*inner_radius, x1=inner_radius, y1=inner_radius, line_color="black")
  fig.add_shape(type="circle", xref="x", yref="y", x0=-1*outer_radius, y0=-1*outer_radius, x1=outer_radius, y1=outer_radius, line_color="black")
  fig.show()

def view_shower_3D(df, eventNum):
  """
  Creates a displays a 3D visualization of the BCal shower in the provided row of the provided dataframe.

  Parameters:
  -----------
  df \: :external:class:`pandas.DataFrame`
      The dataframe containing events in the GlueX BCal. This dataframe should have 7 columns: "module","sector", "layer", "pulse_integral_up", "pulse_integral_down", "t_up", and "t_down". Each of these columns should contain a list, and all of these lists should be the same length across columns. Each element of each list corresponds to a single hit in the BCal.
  eventNum \: int
      The row of the dataframe for which you want to plot the shower.

  """
  def closest_hit(prev_layer,x,y,z):
    prev_distance = 10000000
    solution = (0,0,0)
    for hit in prev_layer.iloc:
      distance = ((x-hit['x'])**2 + (y-hit['y'])**2 + (z-hit['z'])**2)**0.5
      if distance < prev_distance:
        solution = (hit['x'],hit['y'],hit['z'])
        prev_distance = distance
    if prev_distance > ((2*x)**2 + (2*y)**2)**0.5:
      return (x,y,z)
    else:
      return solution
  event = df.iloc[eventNum]
  inner_radius = 65
  inner_module_width = 8.52
  fig = go.FigureWidget()
  fig.update_layout(xaxis_range=[-100,100], yaxis_range=[-100,100],width=520,height=500,xaxis_title="X", yaxis_title="Y")
  points = {'x':[],'y':[],'z':[]}
  energies = []
  for i in range(len(event['sector'])):
    angle = event['module'][i] * -1*np.pi/24 + np.pi
    unit_radius_vector = np.array([np.cos(angle),np.sin(angle)])
    unit_vector = np.array([np.cos(angle+np.pi/2),np.sin(angle+np.pi/2)])
    inner_radius_vector = inner_radius * np.array([np.cos(angle),np.sin(angle)])
    origin = inner_radius_vector - inner_module_width/2 * unit_vector
    position = origin + unit_vector * (8.51*(event['sector'][i]-0.5)/4 + 1.475*(event['sector'][i]-0.5)/2*((event['layer'][i]+0.5)**2-event['layer'][i]-0.5)/20)
    position += 22.46*((event['layer'][i]+0.5)**2-event['layer'][i]-0.5)/20 *unit_radius_vector
    points['x'].append(position[0])
    points['y'].append(position[1])
    z = 1/30*16.75*(event['t_up'][i]-event['t_down'][i])
    points['z'].append(z)
    energies.append(2/130000*(event['pulse_integral_up'][i]*np.exp(z/525)+event['pulse_integral_down'][i]*np.exp((390-z)/525)))
  each_hit_df = pd.DataFrame()
  each_hit_df['module'] = event['module']
  each_hit_df['sector'] = event['sector']
  each_hit_df['layer'] = event['layer']
  each_hit_df['x'] = points['x']
  each_hit_df['y'] = points['y']
  each_hit_df['z'] = points['z']
  each_hit_df['E'] = energies
  each_hit_df['opacity'] = 0.5*(np.array(energies)-min(energies))/(max(energies)-min(energies))+0.5
  previous_layer_hits = pd.DataFrame()
  for i in range(1,5):
    layer_hits = each_hit_df.loc[each_hit_df['layer']==i]
    if i != 1 and len(previous_layer_hits.index)>0 and len(layer_hits.index)>0:
      for hit in layer_hits.iloc:
        angle = hit['module'] * -1*np.pi/24 + np.pi
        unit_radius_vector = np.array([np.cos(angle),np.sin(angle)])
        unit_vector = np.array([np.cos(angle+np.pi/2),np.sin(angle+np.pi/2)])
        inner_radius_vector = inner_radius * np.array([np.cos(angle),np.sin(angle)])
        origin = inner_radius_vector - inner_module_width/2 * unit_vector
        position = origin + unit_vector * (8.51*(hit['sector']-0.5)/4 + 1.475*(hit['sector']-0.5)/2*((hit['layer']+0.5)**2-hit['layer']-0.5)/20)
        position += 22.46*((hit['layer']+0.5)**2-hit['layer']-0.5)/20 *unit_radius_vector
        x2 = position[0]
        y2= position[1]
        z2 = hit['z']
        x1,y1,z1 = closest_hit(each_hit_df.loc[each_hit_df['layer']<i],x2,y2,z2)
        fig.add_scatter3d(x=[x1,x2],y=[y1,y2],z=[z1,z2],mode='lines',line={'color':"rgba(99,110,250,"+str(round(hit['opacity'],2))+")",'width':8})
    previous_layer_hits = layer_hits
  fig.update_layout(showlegend=False)
  fig.show()