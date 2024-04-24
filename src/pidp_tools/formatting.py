def format_labels(particle):
  """
  Formats a particle label from a plain string to a LaTeX string.

  Parameters
  ----------
  particle \: str
      The particle label to be converted into a LaTeX string.

  Returns
  -------
  str
      The LaTeX formatted label, or the original label if not a valid particle type.
  """
  if particle == "Photon":
    return r"$\gamma$"
  elif particle == "KLong":
    return r"$K_L$"
  elif particle == "Neutron":
    return r"$n$"
  elif particle == "Proton":
    return r"$p$"
  elif particle == "K+":
    return r"$K+$"
  elif particle == "Pi+":
    return r"$\pi+$"
  elif particle == "AntiMuon":
    return r"$\overline{\mu}$"
  elif particle == "Positron":
    return r"$\overline{e}$"
  elif particle == "AntiProton":
    return r"$\overline{p}$"
  elif particle == "K-":
    return r"$K-$"
  elif particle == "Pi-":
    return r"$\pi-$"
  elif particle == "Muon":
    return r"$\mu$"
  elif particle == "Electron":
    return r"$e$"
  else:
    return particle
  
def from_pdg(ptype):
  """
  Finds the string label of the particle based on a Particle Data Group (PDG) number.

  Parameters
  ----------
  ptype \: int
      A PDG number corresponding to a particle in the following list\: ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]

  Returns
  -------
  str
      A string of the name of the particle that corresponds to the provided pdg number.
  """
  if int(ptype) == 11:
    return "Electron"
  elif int(ptype) == -11:
    return "Positron"
  elif int(ptype) == 13:
    return "Muon"
  elif int(ptype) == -13:
    return "AntiMuon"
  elif int(ptype) == 211:
    return "Pi+"
  elif int(ptype) == -211:
    return "Pi-"
  elif int(ptype) == 321:
    return "K+"
  elif int(ptype) == -321:
    return "K-"
  elif int(ptype) == 2212:
    return "Proton"
  elif int(ptype) == -2212:
    return "AntiProton"
  elif int(ptype) == 2112:
    return "Neutron"
  elif int(ptype) == 22:
    return "Photon"
  elif int(ptype) == 130:
    return "KLong"
