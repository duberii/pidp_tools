import plotly.graph_objects as go
class vector():
  def __init__(self,x_val, y_val, z_val):
    x = [0, x_val]
    y = [0, y_val]
    z = [0, z_val]
    data = [go.Scatter3d(x=x, y=y, z=z,mode="lines",type="scatter3d",hoverinfo="none",line=dict(color="blue",width=3),name="Vector"),go.Cone(x=[x[1]],y=[y[1]],z=[z[1]],u=[0.3 * (x[1] - x[0])],v=[0.3 * (y[1] - y[0])],w=[0.3 * (z[1] - z[0])],anchor="tip",hoverinfo="none",colorscale=[[0, "blue"], [1, "blue"]],showscale=False),go.Scatter3d(x=[0],y=[0],z=[0],mode="markers",marker=dict(size=5,color="red"),name="Origin")]
    layout = go.Layout(scene=dict(aspectmode="cube",xaxis=dict(range=[-10,10]),yaxis=dict(range=[-10,10]),zaxis=dict(range=[-10,10])))
    self.fig = go.Figure(data=data, layout=layout)
  def add_vector(self, x_val, y_val, z_val):
    x = [0, x_val]
    y = [0, y_val]
    z = [0, z_val]
    self.fig.add_trace(go.Scatter3d(x=x, y=y, z=z,mode="lines",type="scatter3d",hoverinfo="none",line=dict(color="blue",width=3),name="Vector"))
    self.fig.add_trace(go.Cone(x=[x[1]],y=[y[1]],z=[z[1]],u=[0.3 * (x[1] - x[0])],v=[0.3 * (y[1] - y[0])],w=[0.3 * (z[1] - z[0])],anchor="tip",hoverinfo="none",colorscale=[[0, "blue"], [1, "blue"]],showscale=False))