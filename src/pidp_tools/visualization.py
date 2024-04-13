import plotly.graph_objects as go
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from ipywidgets import interact, interactive, fixed, interact_manual, widgets, HBox, VBox, Layout
import ipywidgets as widgets

def plot_vector(x_val,y_val,z_val):
  x = [0, x_val]
  y = [0, y_val]
  z = [0, z_val]
  data = [go.Scatter3d(x=x, y=y, z=z,mode="lines",type="scatter3d",hoverinfo="none",line=dict(color="blue",width=3),name="Vector"),go.Cone(x=[x[1]],y=[y[1]],z=[z[1]],u=[0.3 * (x[1] - x[0])],v=[0.3 * (y[1] - y[0])],w=[0.3 * (z[1] - z[0])],anchor="tip",hoverinfo="none",colorscale=[[0, "blue"], [1, "blue"]],showscale=False),go.Scatter3d(x=[0],y=[0],z=[0],mode="markers",marker=dict(size=5,color="red"),name="Origin")]
  layout = go.Layout(scene=dict(aspectmode="cube",xaxis=dict(range=[-10,10]),yaxis=dict(range=[-10,10]),zaxis=dict(range=[-10,10])))
  fig = go.Figure(data=data, layout=layout)
  return fig.show()

def x_coord(r_curvature, theta_curvature,center_x, t,q):
  return r_curvature*np.cos(theta_curvature-q*t) + center_x

def y_coord(r_curvature,theta_curvature,center_y, t,q):
  return r_curvature*np.sin(theta_curvature-q*t)+ center_y

def degrees_to_radians(degrees):
  return degrees*np.pi/180

def get_wire_positions():
  !wget -q -O CentralDC_HDDS.xml https://github.com/JeffersonLab/hdds/raw/master/CentralDC_HDDS.xml
  tree =ET.parse("CentralDC_HDDS.xml")
  root = tree.getroot()
  n_wires_per_ring = [0]
  wire_to_ring = [0]
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
  wire_positions_df = pd.DataFrame.from_records(ar)
  return wire_positions_df, positionsMatrix


def ring_wire_to_x_y(rings, wires,s, dataframe=False):
  if dataframe:
    df = pd.DataFrame()
    df['ring'] = rings
    df['wire'] = wires
    if not isinstance(s,float) and not isinstance(s, int):
      df['x'] = [positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s[i]) + s[i]* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))]
      df['y'] = [positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s[i]) + s[i]* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))]
    else:
      df['x'] = [positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s) + s* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))]
      df['y'] = [positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s) + s* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))]
    return df
  else:
    if not isinstance(s,float) and not isinstance(s, int):
      new_x = np.array([positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s[i]) + s[i]* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))])
      new_y = np.array([positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s[i]) + s[i]* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))])
    else:
      new_x = np.array([positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][0]*(1-s) + s* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][0] for i in range(len(wires))])
      new_y = np.array([positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][0][1]*(1-s) + s* positionsMatrix[int(rings[i]-1)][int(wires[i]-1)][1][1] for i in range(len(wires))])
    return new_x, new_y

def get_charge(ptype):
  if ptype in ["Electron", "Muon", "Pi-", "K-",'AntiProton']:
    return -1
  elif ptype in ["Positron","AntiMuon", "Pi+", "K+", "Proton"]:
    return 1
  else:
    return 0

class CDC_plot():
  def __init__(self, title, event=0, showTrack = False, showHits = False, showlegend=False):
    self.figure = go.FigureWidget()
    self.figure.update_layout(xaxis_range=[-60,60], yaxis_range=[-60,60],width=500,height=500,showlegend=showlegend,title=title,xaxis_title="X", yaxis_title="Y")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-10.5, y0=-10.5, x1=10.5, y1=10.5, line_color="black")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-55, y0=-55, x1=55, y1=55, line_color="black")
    if not isinstance(event, int):
      showHits = True
      self.charge = get_charge(event["particle"])
      self.px = event['px']
      self.py = event['py']
      self.pz = event['pz']
      self.vz = event['vz']
      self.ring = event['ring']
      self.straw = event['straw']
      df = ring_wire_to_x_y(self.ring,self.straw,75,dataframe=True)
      self.figure.add_scatter(x=df['x'].to_numpy(),y=df['y'].to_numpy(),mode='markers',text = ["Ring: " + str(row['ring']) + "\n Wire: " + str(row['wire']) for row in df.iloc],hoverinfo='text',marker={"size":3})
    if showTrack:
      self.figure.add_scatter(mode='lines',name='track',marker={"color":'red'})
  def show(self):
    self.figure.show()
  def update_layout(self, *args, **kwargs):
    self.figure.update_layout(*args, **kwargs)
  def update(self, *args, **kwargs):
    self.figure.update(*args, **kwargs)
  def display(self,*args,**kwargs):
    widget=interactive(*args,**kwargs)
    controls = HBox(widget.children[:-1], layout = Layout(flex_flow='row wrap'))
    output = widget.children[-1]
    display(VBox([controls, output]))
  def update_figure(self, z0=0,px=0,py=0,charge=1,pz=0,t_max = np.pi, show_track = True, show_hits = True):
    B=1.7
    p = (px**2 + py**2 + pz**2)**0.5
    try:
      z_angle = np.arcsin(pz/p)
    except:
      z_angle = 0
    center_x = charge*330*py/B
    center_y = -1*charge*330*px/B
    r_curvature = (center_x**2 + center_y**2)**0.5
    theta_curvature = np.arctan2(-1*center_y,-1*center_x)
    if show_hits:
      new_x,new_y = ring_wire_to_x_y(self.ring,self.straw,z0/175)
      theta_hits= np.arctan2(new_y-center_y,new_x-center_x)
      for i in range(len(theta_hits)):
        theta_hits[i]-= theta_curvature
        theta_hits[i] = theta_hits[i] * -1 * charge
        if theta_hits[i] <0:
          theta_hits[i] += 2*np.pi
      for i in range(len(self.ring)):
        new_x,new_y = ring_wire_to_x_y(self.ring,self.straw,[(z0 +r_curvature*theta_hits[i]*np.tan(z_angle))/175 for i in range(len(self.ring))])
    if show_track:
      t_curve = np.linspace(0,t_max,1000)
      x_track = x_coord(r_curvature,theta_curvature,center_x,t_curve,charge)
      y_track = y_coord(r_curvature,theta_curvature,center_y,t_curve,charge)
    with fig.figure.batch_update():
      if show_track:
        fig.figure.data[len(fig.figure.data)-1]['x']=x_track
        fig.figure.data[len(fig.figure.data)-1]['y']=y_track
      fig.update_layout(width=500,height=500,showlegend=False)
      if show_hits:
        fig.figure.data[0]['x']=new_x
        fig.figure.data[0]['y']=new_y
    fig.show()
  def distance_to_curve(self,x,y):
    center_x = self.charge*330*self.py/1.7
    center_y = -330*self.charge*self.px/1.7
    r_curvature = (center_x**2 + center_y**2)**0.5
    distance_to_nearest_curve_point = abs(r_curvature-((x-center_x)**2+(y-center_y)**2)**0.5)
    return distance_to_nearest_curve_point
  def calculate_rmse(self,x_points, y_points):
    mse = sum([self.distance_to_curve(x_points[j],y_points[j])**2 for j in range(len(x_points))])/len(x_points)
    rmse = mse**0.5
    return rmse
  def calc_z0(self, quiet=False, return_rmse=False):
    rmses = {}
    center_x = self.charge*330*self.py/B
    center_y = -330*self.charge*self.px/B
    for z in np.linspace(0,150,151):
      new_x, new_y = ring_wire_to_x_y(self.ring, self.straw,z/175)
      rmses[z] = self.calculate_rmse(new_x,new_y)
    if quiet and return_rmse:
      return min(rmses, key=rmses.get), min(list(rmses.values()))
    elif quiet:
      return min(rmses, key=rmses.get)
    elif return_rmse:
      min(list(rmses.values()))
    else:
      print("The true z value is about " + str(min(rmses, key=rmses.get)))

  def show_answer(self, charge=True,px=True,py=True,pz=False,z0=False,but="", every=False):
    if z0 or all([every,but != "z0"]):
      self.calc_z0()
    if charge or all([every,but != "charge"]):
      print("The charge of the particle is " + str(self.charge))
    if px or all([every,but != "px"]):
      print("The x component of momentum is " + str(round(self.px, 2)))
    if py or all([every,but != "py"]):
      print("The y component of momentum is " + str(round(self.py, 2)))
    if pz or all([every,but != "pz"]):
      print("The z component of momentum is " + str(round(self.pz, 2)))

class wire_plot():
  def __init__(self,rings,title, showlegend=True):
    self.figure = go.FigureWidget()
    self.figure.update_layout(xaxis_range=[-60,60], yaxis_range=[-60,60],width=500,height=500,showlegend=showlegend,title=title,xaxis_title="X", yaxis_title="Y")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-10.5, y0=-10.5, x1=10.5, y1=10.5, line_color="black")
    self.figure.add_shape(type="circle", xref="x", yref="y", x0=-55, y0=-55, x1=55, y1=55, line_color="black")
    for ring in rings:
      self.figure.add_scatter(mode='markers',name="Ring " + str(ring),marker={"size":3},hoverinfo='text',text = ['Ring: ' + str(ring) + "\n Wire: " + str(k) for k in range(len(positionsMatrix[int(ring)-1]))])
    self.data = []
  def show(self):
    self.figure.show()
  def update_layout(self, *args, **kwargs):
    self.figure.update_layout(*args, **kwargs)
  def update(self, *args, **kwargs):
    self.figure.update(*args, **kwargs)
  def display(self,*args,**kwargs):
    widget=interactive(*args,**kwargs)
    controls = HBox(widget.children[:-1], layout = Layout(flex_flow='row wrap'))
    output = widget.children[-1]
    display(VBox([controls, output]))
  def update_figure(self,z=0, rings=[]):
    s= z/175
    with fig.figure.batch_update():
      for i in range(len(rings)):
        fig.figure.data[i]['x']=[(1-s)*wire[0][0]+ s*wire[1][0] for wire in positionsMatrix[int(rings[i])-1]]
        fig.figure.data[i]['y']=[(1-s)*wire[0][1]+ s*wire[1][1] for wire in positionsMatrix[int(rings[i])-1]]
    fig.show()

class interactive_image():
  def __init__(self,image_path="",x_max = 1, y_max=2*10**-5):
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
    self.figure.show()
  def update_layout(self, *args, **kwargs):
    self.figure.update_layout(*args, **kwargs)
  def update(self, *args, **kwargs):
    self.figure.update(*args, **kwargs)
  def display(self,*args,**kwargs):
    widget=interactive(*args,**kwargs)
    controls = HBox(widget.children[:-1], layout = Layout(flex_flow='row wrap'))
    output = widget.children[-1]
    display(VBox([controls, output]))
  def update_figure(self,a=1, b=1, c=1):
    with self.figure.batch_update():
      self.figure.data[0]['y']=10**(math.floor(math.log(self.y_max,10))-1)*(np.exp(a*np.linspace(0,self.x_max,1000)+b)+c)
    self.show()