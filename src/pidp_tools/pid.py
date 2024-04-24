import numpy as np

def dEdx_CDC(particle_df):
  """
  Identifies particles based on their ionization energy loss in the GlueX Central Drift Chamber (CDC).

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", and "CDC dEdx".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on CDC ionization energy loss cuts.
  """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (particle['CDC dEdx'] *10**6> np.exp(-2.5*p+2.3)+1.4)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['CDC dEdx'] *10**6> np.exp(-4.9*p+2.3)+1.45) & (particle['CDC dEdx'] *10**6< np.exp(-4.36*p+2.3)+2.3)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['CDC dEdx'] *10**6> 1.59) & (particle['CDC dEdx'] *10**6< 2.49)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['CDC dEdx'] *10**6> 2.25) & (particle['CDC dEdx'] *10**6< 3.14)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions
    

def dEdx_FDC(particle_df):
  """
  Identifies particles based on their ionization energy loss in the GlueX Forward Drift Chambers (FDC).

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", and "FDC dEdx".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on FDC ionization energy loss cuts.
  """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (particle['FDC dEdx'] *10**6> np.exp(-3.55*p+3.7)+1.3)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['FDC dEdx'] *10**6> np.exp(-4.4*p+2)+1.3) & (particle['FDC dEdx'] *10**6< np.exp(-4.4*p+2)+2.3)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['FDC dEdx'] *10**6> 1.45) & (particle['FDC dEdx'] *10**6 < 2.6)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['FDC dEdx'] *10**6> 2) & (particle['FDC dEdx'] *10**6 < 3)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions
  
def dEdx_SC(particle_df):
  """
  Identifies particles based on their ionization energy loss in the GlueX Start Counter (SC).

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", and "SC dEdx".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on SC ionization energy loss cuts.
    """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (particle['SC dEdx'] *10**3> np.exp(-4*p+3.2)+1.6)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['SC dEdx'] *10**3> np.exp(-4.4*p+2.7)+1.2) & (particle['SC dEdx'] *10**3< np.exp(-4.6*p+3.1)+2.2)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['SC dEdx'] *10**3> 1.3) & (particle['SC dEdx'] *10**3< 2.3)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['SC dEdx'] *10**3> 1.3) & (particle['SC dEdx'] *10**3< 2.3)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions
  
def dEdx_TOF(particle_df):
  """
  Identifies particles based on their ionization energy loss in the GlueX Time of Flight (TOF) detector.

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", and "TOF dEdx".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on TOF ionization energy loss cuts.
  """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (particle['TOF dEdx'] *10**4> np.exp(-2.2*p+3.9)+8.2)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['TOF dEdx'] *10**4> np.exp(-0.4*p+0.8)+7.6) & (particle['TOF dEdx'] *10**4< np.exp(-0.4*p+1.4)+10)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['TOF dEdx'] *10**4> 8.6) & (particle['TOF dEdx'] *10**4< 11)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['TOF dEdx'] *10**4> 8.1) & (particle['TOF dEdx'] *10**4< 11)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions
  
def TOF_TOF(particle_df):
  """
  Identifies particles based on their time of flight, as measured by the GlueX Time of Flight (TOF) detector.

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", "TOF Time of Flight", and "TOF Calculated Time of Flight".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on TOF TOF cuts.
  """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  particle['TOF TOF Diff'] = np.abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight'])
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (particle['TOF TOF Diff'] <= 0.6)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['TOF TOF Diff'] <= 0.3)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['TOF TOF Diff'] <= 0.5)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['TOF TOF Diff'] <= 0.5)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions
  
def TOF_BCal(particle_df):
  """
  Identifies particles based on their time of flight, as measured by the GlueX Barrel Calorimeter (BCal).

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", "tShower", and "BCal Calculated Time of Flight".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on BCal TOF cuts.
    """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  particle['BCal TOF Diff'] = np.abs(particle['tShower']-particle['BCal Calculated Time of Flight'])
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & ( particle['BCal TOF Diff'] <= 1)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['BCal TOF Diff'] <= 0.75)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['BCal TOF Diff'] <= 1)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['BCal TOF Diff'] <= 1)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions

def TOF_SC(particle_df):
  """
  Identifies particles based on their time of flight, as measured by the GlueX Start Counter (SC).

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", "SC Time of Flight", and "SC Calculated Time of Flight".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on SC TOF cuts.
  """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  particle['SC TOF Diff'] = np.abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight'])
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (particle['SC TOF Diff'] <= 2.5)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['SC TOF Diff'] <= 2.5)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['SC TOF Diff'] <= 2.5)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['SC TOF Diff'] <= 2.5)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions

def TOF_FCal(particle_df):
  """
  Identifies particles based on their time of flight, as measured by the GlueX Forward Calorimeter (FCal).

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", "tShower", and "FCal Calculated Time of Flight".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on FCal TOF cuts.
  """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  particle['FCal TOF Diff'] = np.abs(particle['tShower']-particle['FCal Calculated Time of Flight'])
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (particle['FCal TOF Diff'] <= 2)
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (particle['FCal TOF Diff'] <= 2.5)
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (particle['FCal TOF Diff'] <= 2)
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (particle['FCal TOF Diff'] <= 2)
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions

def TOF(particle_df):
  """
  Identifies particles based on their best-available time of flight, as measured by the GlueX Barrel Calorimeter (BCal), Time of Flight (TOF) detector, Forward Calorimeter (FCal), or the Start Counter (SC).

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", "tShower", "BCal Calculated Time of Flight", "TOF Time of Flight", "TOF Calculated Time of Flight", "FCal Calculated Time of Flight", "SC Time of Flight", and "SC Calculated Time of Flight".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on best-available TOF cuts.
    """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  particle['BCal TOF Diff'] = np.abs(particle['tShower']-particle['BCal Calculated Time of Flight'])
  particle['TOF TOF Diff'] = np.abs(particle['TOF Time of Flight']-particle['TOF Calculated Time of Flight'])
  particle['FCal TOF Diff'] = np.abs(particle['tShower']-particle['FCal Calculated Time of Flight'])
  particle['SC TOF Diff'] = np.abs(particle['SC Time of Flight']-particle['SC Calculated Time of Flight'])
  hasBCal = particle['BCal Calculated Time of Flight'].notna()
  hasTOF = particle['TOF Time of Flight'].notna()
  hasFCal = particle['FCal Calculated Time of Flight'].notna()
  hasSC = particle['SC Time of Flight'].notna()
  p = (particle['px']**2 + particle['py']**2 + particle['pz']**2)**0.5
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & ((hasBCal & (particle['BCal TOF Diff'] <= 1)) | (~hasBCal & hasTOF & (particle['TOF TOF Diff'] <= 0.6)) | (~hasBCal & ~hasTOF & hasFCal & (particle['FCal TOF Diff'] <= 2)) | (~hasBCal & ~hasTOF & ~hasFCal & hasSC & (particle['SC TOF Diff'] <= 2.5)))
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & ((hasBCal & (particle['BCal TOF Diff'] <= 0.75)) | (~hasBCal & hasTOF & (particle['TOF TOF Diff'] <= 0.3)) | (~hasBCal & ~hasTOF & hasFCal & (particle['FCal TOF Diff'] <= 2.5)) | (~hasBCal & ~hasTOF & ~hasFCal & hasSC & (particle['SC TOF Diff'] <= 2.5)))
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & ((hasBCal & (particle['BCal TOF Diff'] <= 1)) | (~hasBCal & hasTOF & (particle['TOF TOF Diff'] <= 0.5)) | (~hasBCal & ~hasTOF & hasFCal & (particle['FCal TOF Diff'] <= 2)) | (~hasBCal & ~hasTOF & ~hasFCal & hasSC & (particle['SC TOF Diff'] <= 2.5)))
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & ((hasBCal & (particle['BCal TOF Diff'] <= 1)) | (~hasBCal & hasTOF & (particle['TOF TOF Diff'] <= 0.5)) | (~hasBCal & ~hasTOF & hasFCal & (particle['FCal TOF Diff'] <= 2)) | (~hasBCal & ~hasTOF & ~hasFCal & hasSC & (particle['SC TOF Diff'] <= 2.5)))
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions
  
def DIRC(particle_df):
  """
  Identifies particles based on their likelihood values produced by the GlueX Detection of Internally Reflected Cherenkov Radiation (DIRC) detector.

  Parameters
  ----------
  particle_df \: :external:class:`pandas.DataFrame`
      A dataframe containing particles to be identified. The supplied dataframe must have the following columns: "px","py","pz", "Hypothesis", "lp", "lk", "lpi", and "lele".

  Returns
  -------
  predictions \: :external:class:`pandas.Series`
      The predictions based on DIRC likelihood cuts.
  """
  particle = particle_df.copy()
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if not isinstance(particle['Hypothesis'][0], str):
    particle['Hypothesis'].apply(lambda x: particle_list[int(x)])
  hasDIRC = particle['lp'].notna()
  maxLik = particle[['lp','lk','lpi','lele']].max()
  isProton = ((particle['Hypothesis'] == "Proton") | (particle['Hypothesis'] == "AntiProton")) & (maxLik == particle['lp']) & hasDIRC
  isK = ((particle['Hypothesis'] == "K+") | (particle['Hypothesis'] == "K-")) & (maxLik == particle['lk']) & hasDIRC
  isPi = ((particle['Hypothesis'] == "Pi+") | (particle['Hypothesis'] == "Pi-")) & (maxLik == particle['lpi']) & hasDIRC
  isElectron = ((particle['Hypothesis'] == "Electron") | (particle['Hypothesis'] == "Positron")) & (maxLik == particle['lele']) & hasDIRC
  particle['matches hypothesis'] = isProton | isK | isPi | isElectron
  def temp_match(row):
    if row['matches hypothesis']:
      return row['Hypothesis']
    else:
      return "No ID"
  predictions = particle.apply(temp_match, axis=1)
  return predictions