import plotly.graph_objects as go
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from ipywidgets import interact, interactive, fixed, interact_manual, widgets, HBox, VBox, Layout
import ipywidgets as widgets
import subprocess
from pidp_tools.analysis import get_charge

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

class wirePositions():
  def __init__(self):
    subprocess.run("wget -q -O CentralDC_HDDS.xml https://github.com/JeffersonLab/hdds/raw/master/CentralDC_HDDS.xml")
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
    self.wire_positions_df = pd.DataFrame.from_records(ar)
    self.positionsMatrix = positionsMatrix
  def position(self, rings, wires,s, dataframe=False):
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

class CDC_plot():
  def __init__(self, title, wire_positions, event=0, showTrack = False, showHits = False, showlegend=False):
    self.wirePositions = wire_positions
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
      df = self.wirePositions.position(self.ring,self.straw,75,dataframe=True)
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
    center_x = self.charge*330*self.py/1.7
    center_y = -330*self.charge*self.px/1.7
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
    with self.figure.batch_update():
      for i in range(len(rings)):
        self.figure.data[i]['x']=[(1-s)*wire[0][0]+ s*wire[1][0] for wire in positionsMatrix[int(rings[i])-1]]
        self.figure.data[i]['y']=[(1-s)*wire[0][1]+ s*wire[1][1] for wire in positionsMatrix[int(rings[i])-1]]
    self.show()

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


def BCal_plot_2D():
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
  def closest_hit(prev_layer,x,y):
    prev_distance = 10000000
    solution = (0,0)
    for hit in prev_layer.iloc:
      distance = ((x-hit['x'])**2 + (y-hit['y'])**2)**0.5
      if distance < prev_distance:
        solution = (hit['x'],hit['y'],hit['opacity'])
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