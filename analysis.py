class interactive_image():
  def __init__(self,image_path=""):
    self.figure = go.FigureWidget()
    self.figure.add_scatter(x=np.linspace(0,1,1000), y=np.zeros(1000),mode='lines',name='Cut',marker={"color":'red'})
    if len(image_path)>0:
      self.figure.add_layout_image({"source":Image.open(image_path),"xref":"x", "yref":"y","x": -0.125,"y":2.25*10**-5, "sizex": 1.25,"sizey":2.5*10**-5,"sizing":"stretch","opacity":1,"layer":"below"})
    self.figure.update_xaxes(range=[0, 1],showgrid=False)
    self.figure.update_yaxes(range=[0, 2*10**-5],showgrid=False)
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
      self.figure.data[0]['y']=10**-6*(np.exp(a*np.linspace(0,1,1000)+b)+c)
    self.show()