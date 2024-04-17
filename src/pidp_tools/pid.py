import numpy as np

def dEdx_CDC_Proton(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['CDC dEdx'] *10**6> np.exp(-2.5*p+2.3)+1.4:
    return True
  else:
    return False

def dEdx_CDC_Kaon(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['CDC dEdx'] *10**6> np.exp(-4.9*p+2.3)+1.45 and particle['CDC dEdx'] *10**6< np.exp(-4.36*p+2.3)+2.3:
    return True
  else:
    return False

def dEdx_CDC_Pion(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['CDC dEdx'] *10**6> 1.59 and particle['CDC dEdx'] *10**6< 2.49:
    return True
  else:
    return False

def dEdx_CDC_Electron(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['CDC dEdx'] *10**6> 2.25 and particle['CDC dEdx'] *10**6< 3.14:
    return True
  else:
    return False

def dEdx_CDC(particle):
  if particle['Hypothesis'] == "Proton" and dEdx_CDC_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and dEdx_CDC_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and dEdx_CDC_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and dEdx_CDC_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and dEdx_CDC_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and dEdx_CDC_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and dEdx_CDC_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and dEdx_CDC_Electron(particle):
    return "Positron"
  else:
    return "No ID"
    

def dEdx_FDC_Proton(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['FDC dEdx'] *10**6> np.exp(-3.55*p+3.7)+1.3:
    return True
  else:
    return False

def dEdx_FDC_Kaon(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['FDC dEdx'] *10**6> np.exp(-4.4*p+2)+1.3 and particle['FDC dEdx'] *10**6< np.exp(-4.4*p+2)+2.3:
    return True
  else:
    return False

def dEdx_FDC_Pion(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['FDC dEdx'] *10**6> 1.45 and particle['FDC dEdx'] *10**6 < 2.6:
    return True
  else:
    return False

def dEdx_FDC_Electron(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['FDC dEdx'] *10**6> 2 and particle['FDC dEdx'] *10**6 < 3:
    return True
  else:
    return False

def dEdx_FDC(particle):
  if particle['Hypothesis'] == "Proton" and dEdx_FDC_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and dEdx_FDC_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and dEdx_FDC_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and dEdx_FDC_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and dEdx_FDC_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and dEdx_FDC_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and dEdx_FDC_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and dEdx_FDC_Electron(particle):
    return "Positron"
  else:
    return "No ID"
  
def SC_dEdx_Proton(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['SC dEdx'] *10**3> np.exp(-4*p+3.2)+1.6:
    return True
  else:
    return False

def SC_dEdx_Kaon(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['SC dEdx'] *10**3> np.exp(-4.4*p+2.7)+1.2 and particle['SC dEdx'] *10**3< np.exp(-4.6*p+3.1)+2.2:
    return True
  else:
    return False

def SC_dEdx_Pion(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['SC dEdx'] *10**3> 1.3 and particle['SC dEdx'] *10**3< 2.3:
    return True
  else:
    return False

def SC_dEdx_Electron(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['SC dEdx'] *10**3> 1.3 and particle['SC dEdx'] *10**3< 2.3:
    return True
  else:
    return False

def dEdx_SC(particle):
  if particle['Hypothesis'] == "Proton" and SC_dEdx_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and SC_dEdx_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and SC_dEdx_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and SC_dEdx_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and SC_dEdx_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and SC_dEdx_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and SC_dEdx_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and SC_dEdx_Electron(particle):
    return "Positron"
  else:
    return "No ID"
  
def TOF_dEdx_Proton(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['TOF dEdx'] *10**4> np.exp(-2.2*p+3.9)+8.2:
    return True
  else:
    return False

def TOF_dEdx_Kaon(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['TOF dEdx'] *10**4> np.exp(-0.4*p+0.8)+7.6 and particle['TOF dEdx'] *10**4< np.exp(-0.4*p+1.4)+10:
    return True
  else:
    return False

def TOF_dEdx_Pion(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['TOF dEdx'] *10**4> 8.6 and particle['TOF dEdx'] *10**4< 11:
    return True
  else:
    return False

def TOF_dEdx_Electron(particle):
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  if particle['TOF dEdx'] *10**4> 8.1 and particle['TOF dEdx'] *10**4< 11:
    return True
  else:
    return False

def dEdx_TOF(particle):
  if particle['Hypothesis'] == "Proton" and TOF_dEdx_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and TOF_dEdx_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and TOF_dEdx_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and TOF_dEdx_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and TOF_dEdx_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and TOF_dEdx_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and TOF_dEdx_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and TOF_dEdx_Electron(particle):
    return "Positron"
  else:
    return "No ID"
  
def TOF_TOF_Proton(particle):
  if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.6:
    return True
  else:
    return False

def TOF_TOF_Kaon(particle):
  if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.3:
    return True
  else:
    return False

def TOF_TOF_Pion(particle):
  if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.5:
    return True
  else:
    return False

def TOF_TOF_Electron(particle):
  if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.5:
    return True
  else:
    return False

def TOF_TOF(particle):
  if particle['Hypothesis'] == "Proton" and TOF_TOF_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and TOF_TOF_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and TOF_TOF_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and TOF_TOF_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and TOF_TOF_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and TOF_TOF_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and TOF_TOF_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and TOF_TOF_Electron(particle):
    return "Positron"
  else:
    return "No ID"
  
def BCal_TOF_Proton(particle):
  if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 1:
    return True
  else:
    return False

def BCal_TOF_Kaon(particle):
  if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 0.75:
    return True
  else:
    return False

def BCal_TOF_Pion(particle):
  if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 1:
    return True
  else:
    return False

def BCal_TOF_Electron(particle):
  if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 1:
    return True
  else:
    return False

def TOF_BCal(particle):
  if particle['Hypothesis'] == "Proton" and BCal_TOF_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and BCal_TOF_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and BCal_TOF_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and BCal_TOF_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and BCal_TOF_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and BCal_TOF_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and BCal_TOF_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and BCal_TOF_Electron(particle):
    return "Positron"
  else:
    return "No ID"

def SC_TOF_Proton(particle):
  if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
    return True
  else:
    return False

def SC_TOF_Kaon(particle):
  if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
    return True
  else:
    return False

def SC_TOF_Pion(particle):
  if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
    return True
  else:
    return False

def SC_TOF_Electron(particle):
  if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
    return True
  else:
    return False

def TOF_SC(particle):
  if particle['Hypothesis'] == "Proton" and SC_TOF_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and SC_TOF_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and SC_TOF_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and SC_TOF_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and SC_TOF_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and SC_TOF_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and SC_TOF_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and SC_TOF_Electron(particle):
    return "Positron"
  else:
    return "No ID"

def FCal_TOF_Proton(particle):
  if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2:
    return True
  else:
    return False

def FCal_TOF_Kaon(particle):
  if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2.5:
    return True
  else:
    return False

def FCal_TOF_Pion(particle):
  if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2:
    return True
  else:
    return False

def FCal_TOF_Electron(particle):
  if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2:
    return True
  else:
    return False

def TOF_FCal(particle):
  if particle['Hypothesis'] == "Proton" and FCal_TOF_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and FCal_TOF_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and FCal_TOF_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and FCal_TOF_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and FCal_TOF_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and FCal_TOF_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and FCal_TOF_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and FCal_TOF_Electron(particle):
    return "Positron"
  else:
    return "No ID"
  
def TOF_Proton(particle):
  if not np.isnan(particle['BCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 1:
      return True
    else:
      return False
  elif not np.isnan(particle['TOF Calculated Time of Flight']):
    if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.6:
      return True
    else:
      return False
  elif not np.isnan(particle['FCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2:
      return True
    else:
      return False
  elif not np.isnan(particle['SC Calculated Time of Flight']):
    if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
      return True
    else:
      return False
  else:
    return False

def TOF_Kaon(particle):
  if not np.isnan(particle['BCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 0.75:
      return True
    else:
      return False
  elif not np.isnan(particle['TOF Calculated Time of Flight']):
    if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.3:
      return True
    else:
      return False
  elif not np.isnan(particle['FCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2.5:
      return True
    else:
      return False
  elif not np.isnan(particle['SC Calculated Time of Flight']):
    if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
      return True
    else:
      return False
  else:
    return False

def TOF_Pion(particle):
  if not np.isnan(particle['BCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 1:
      return True
    else:
      return False
  elif not np.isnan(particle['TOF Calculated Time of Flight']):
    if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.5:
      return True
    else:
      return False
  elif not np.isnan(particle['FCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2:
      return True
    else:
      return False
  elif not np.isnan(particle['SC Calculated Time of Flight']):
    if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
      return True
    else:
      return False
  else:
    return False

def TOF_Electron(particle):
  if not np.isnan(particle['BCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['BCal Calculated Time of Flight']) <= 1:
      return True
    else:
      return False
  elif not np.isnan(particle['TOF Calculated Time of Flight']):
    if abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight']) <= 0.5:
      return True
    else:
      return False
  elif not np.isnan(particle['FCal Calculated Time of Flight']):
    if abs(particle['tShower']-particle['FCal Calculated Time of Flight']) <= 2:
      return True
    else:
      return False
  elif not np.isnan(particle['SC Calculated Time of Flight']):
    if abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight']) <= 2.5:
      return True
    else:
      return False
  else:
    return False

def Combined_TOF(particle):
  if particle['Hypothesis'] == "Proton" and TOF_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and TOF_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and TOF_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and TOF_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and TOF_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and TOF_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and TOF_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and TOF_Electron(particle):
    return "Positron"
  else:
    return "No ID"
  
def DIRC_Proton(particle):
  if particle['lp'] > max(particle['lk'], particle['lpi'],particle['lele']):
    return True
  else:
    return False

def DIRC_Kaon(particle):
  if particle['lk'] > max(particle['lp'], particle['lpi'],particle['lele']):
    return True
  else:
    return False

def DIRC_Pion(particle):
  if particle['lpi'] > max(particle['lk'], particle['lp'],particle['lele']):
    return True
  else:
    return False

def DIRC_Electron(particle):
  if particle['lele'] > max(particle['lk'], particle['lp'],particle['lpi']):
    return True
  else:
    return False

def DIRC_DIRC(particle):
  if particle['Hypothesis'] == "Proton" and DIRC_Proton(particle):
    return "Proton"
  elif particle['Hypothesis'] == "AntiProton" and DIRC_Proton(particle):
    return "AntiProton"
  elif particle['Hypothesis'] == "K+" and DIRC_Kaon(particle):
    return "K+"
  elif particle['Hypothesis'] == "K-" and DIRC_Kaon(particle):
    return "K-"
  elif particle['Hypothesis'] == "Pi+" and DIRC_Pion(particle):
    return "Pi+"
  elif particle['Hypothesis'] == "Pi-" and DIRC_Pion(particle):
    return "Pi-"
  elif particle['Hypothesis'] == "Electron" and DIRC_Electron(particle):
    return "Electron"
  elif particle['Hypothesis'] == "Positron" and DIRC_Electron(particle):
    return "Positron"
  else:
    return "No ID"