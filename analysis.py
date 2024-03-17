def install_ROOT():
  !wget -q https://github.com/MohamedElashri/ROOT/releases/download/ubuntu/root_v6.28.04_Ubuntu_20.04.zip
  !unzip -o -qq /content/root_v6.28.04_Ubuntu_20.04.zip
  !apt-get -qq install git dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev tar gfortran subversion &> /dev/null 2>&1
  !rm -f root_v6.28.04_Ubuntu_20.04.zip
  !wget -q http://archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2_amd64.deb
  !sudo dpkg -i libssl1.1_1.1.1f-1ubuntu2_amd64.deb &> /dev/null 2>&1
  !rm -f libssl1.1_1.1.1f-1ubuntu2_amd64.deb
  import sys
  sys.path.append("/content/root_build/")
  sys.path.append("/content/root_build/bin/")
  sys.path.append("/content/root_build/include/")
  sys.path.append("/content/root_build/lib/")
  import ctypes
  ctypes.cdll.LoadLibrary('/content/root_build/lib//libCore.so')
  ctypes.cdll.LoadLibrary('/content/root_build/lib//libThread.so')
  ctypes.cdll.LoadLibrary('/content/root_build/lib//libTreePlayer.so')

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