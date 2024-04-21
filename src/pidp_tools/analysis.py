import numpy as np
import matplotlib.pyplot as plt
from pidp_tools.formatting import format_labels
import sklearn.metrics
import joblib
import tqdm
import pandas as pd
import warnings

def install_ROOT():
  import subprocess
  try:
    import ROOT
  except:
    try:
      subprocess.run(["wget","-q","https://github.com/MohamedElashri/ROOT/releases/download/ubuntu/root_v6.28.04_Ubuntu_20.04.zip"])
      subprocess.run(["unzip", "-o", "-qq", "/content/root_v6.28.04_Ubuntu_20.04.zip"])
      subprocess.run(["apt-get", "-qq", "install", "git", "dpkg-dev", "cmake", "g++", "gcc", "binutils", "libx11-dev", "libxpm-dev", "libxft-dev", "libxext-dev", "tar", "gfortran", "subversion", "&>", "/dev/null", "2>&1"])
      subprocess.run(["rm", "-f", "root_v6.28.04_Ubuntu_20.04.zip"])
      subprocess.run(["wget", "-q", "http://archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2_amd64.deb"])
      subprocess.run(["sudo", "dpkg", "-i", "libssl1.1_1.1.1f-1ubuntu2_amd64.deb"])
      subprocess.run(["rm", "-f", "libssl1.1_1.1.1f-1ubuntu2_amd64.deb"])
      import sys
      sys.path.append("/content/root_build/")
      sys.path.append("/content/root_build/bin/")
      sys.path.append("/content/root_build/include/")
      sys.path.append("/content/root_build/lib/")
      import ctypes
      ctypes.cdll.LoadLibrary('/content/root_build/lib//libCore.so')
      ctypes.cdll.LoadLibrary('/content/root_build/lib//libThread.so')
      ctypes.cdll.LoadLibrary('/content/root_build/lib//libTreePlayer.so')
    except:
      raise OSError("Unable to install ROOT. This install was designed specifically for use in google colab. Installing on personal computers is discouraged due to the file size.")

def get_charge(ptype):
  if ptype in ["Electron", "Muon", "Pi-", "K-",'AntiProton',3,4,5,6,7]:
    return -1
  elif ptype in ["Positron","AntiMuon", "Pi+", "K+", "Proton",8,9,10,11,12]:
    return 1
  else:
    return 0

def round_accuracies(num):
  new_num = round(num,2)
  if new_num == 0.00:
    return 0
  else:
    return new_num

class ConfusionMatrix():
  def __init__(self, labels, predictions, target='Generated As',title="", purity=False, label_selection="necessary"):
    self.target = 'Generated As'
    self.title= title
    self.purity = purity
    self.label_selection = label_selection

    particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
    particle_array = np.array(particle_list)

    if isinstance(labels[0],str):
      labels = [particle_list.index(label) for label in labels]
    if isinstance(predictions[0],str):
      predictions = [particle_list.index(prediction) for prediction in predictions]

    labels = [int(i) for i in labels]
    predictions = [int(i) for i in predictions]
    self.calculate_matrix(labels, predictions)
    self.display_matrix(title)
  
  @classmethod
  def from_estimator(cls, estimator, df, target='Generated As', title="", purity=False, label_selection="necessary"):

    # Initialize variables
    particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
    dataset = df.copy().reset_index(drop=True)
    predictions = []
    identities = []

    #Ensure hypothesis and generated as columns are strings for use in PID functions

    if isinstance(df['Hypothesis'][0],np.int64) or isinstance(df['Hypothesis'][0],int):
      dataset['Hypothesis']=dataset['Hypothesis'].apply(lambda x: particle_list[x])
    if isinstance(df['Generated As'][0],np.int64) or isinstance(df['Generated As'][0],int):
      dataset['Generated As']=dataset['Generated As'].apply(lambda x: particle_list[x])

    try:
      predictions_full = estimator(dataset)
    except:
      predictions_full = dataset.apply(estimator,axis=1)

    #Converts predictions, as well as the hypothesis and generated as columns, back to integers.

    predictions_full = predictions_full.apply(particle_list.index)
    dataset['Hypothesis'] = dataset['Hypothesis'].apply(particle_list.index)
    dataset['Generated As'] = dataset['Generated As'].apply(particle_list.index)

    #Analyzes the predictions using the hypothesis scheme

    nRows = len(dataset.index)
    starting_index = 0
    dataset['is matched'] = dataset['Hypothesis'] == predictions_full
    index_list = []
    event_nos = []
    i = 0
    while starting_index < nRows - 1 :
      ending_index = starting_index + dataset['Number of Hypotheses'][starting_index]
      index_list.append((starting_index, ending_index))
      event_nos.extend([i for _ in range(starting_index, ending_index)])
      i+= 1
      starting_index = ending_index
    number_of_events = i
    dataset['eventno'] = event_nos
    reduced_dataset = dataset.loc[dataset['is matched']]
    grouped = reduced_dataset.sample(frac=1)[['Hypothesis','eventno']].groupby('eventno')
    predictions = grouped.head(1).set_index('eventno').sort_index().reindex(list(range(number_of_events)),fill_value=13)['Hypothesis'].to_list()
    identities = [dataset[target][starting_index] for starting_index, ending_index in index_list]

    confusion_matrix = cls(identities, predictions, target=target,title=title, purity=purity, label_selection=label_selection)
    return confusion_matrix
  
  @classmethod
  def from_model(cls, model, df, target="Generated As", title="", purity=False, match_hypothesis=False, label_selection="charge"):

    particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
    dataset = df.copy().reset_index(drop=True)
    predictions = []
    identities = []
    data_to_test = dataset[[column for column in model.feature_names_in_]]

    if isinstance(df['Hypothesis'][0],str):
      dataset['Hypothesis']=dataset['Hypothesis'].apply(particle_list.index)
    if isinstance(df['Generated As'][0],str):
      dataset['Generated As']=dataset['Generated As'].apply(particle_list.index)

    nRows = len(dataset.index)
    starting_index = 0

    index_list = []
    event_nos = []
    i = 0
    while starting_index < nRows - 1 :
      ending_index = starting_index + dataset['Number of Hypotheses'][starting_index]
      index_list.append((starting_index, ending_index))
      event_nos.extend([i for _ in range(starting_index, ending_index)])
      i+= 1
      starting_index = ending_index
    number_of_events = i
    dataset['eventNo'] = event_nos

    if match_hypothesis:
      temp_predictions = model.predict_proba(data_to_test)
      dataset['Confidence'] = [max(probs) for probs in temp_predictions]
      dataset['Prediction'] = np.argmax(temp_predictions, axis=1)

      identities_grouped = dataset[['Generated As','Prediction','eventNo']].groupby('eventNo')
      identities = identities_grouped['Generated As'].head(1).to_list()

      matches_hypotheses_bool_list = (dataset['Prediction'] == dataset['Hypothesis']) | (dataset['Hypothesis'] == 13)
      matching_hypotheses = dataset.loc[matches_hypotheses_bool_list]
      grouped_df = matching_hypotheses[['Generated As','Prediction','Confidence','eventNo']].groupby('eventNo')
      max_confidence_indices = grouped_df['Confidence'].idxmax()
      predictions_temp = dataset['Prediction','eventNo'].iloc[max_confidence_indices]
      
      predictions = predictions_temp.set_index('eventNo').sort_index().reindex(list(range(number_of_events)),fill_value=13).to_list
      identities = grouped_df[target].head(1).to_list()
    else:
      dataset['Prediction'] = model.predict(data_to_test)
      grouped_df = dataset[['Generated As','Prediction','eventNo']].groupby('eventNo')
      identities = grouped_df['Generated As'].head(1).to_list()
      predictions = grouped_df['Prediction'].agg(lambda x: x.value_counts().index[0]).to_list()

    confusion_matrix = cls(identities, predictions, target=target,title=title, purity=purity, label_selection=label_selection)
    return confusion_matrix

  def calculate_matrix(self, labels, predictions):
    temp_confusion_matrix = np.zeros((13,14))
    particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
    particle_array = np.array(particle_list)

    match self.label_selection:
      case 'charge':
        included_particles = list(set(list(labels) + list(predictions)))
        contains_neutral_particles = np.any([i in included_particles for i in [0,1,2]])
        contains_charged_particles = np.any([i in included_particles for i in [3,4,5,6,7,8,9,10,11,12]])
        if contains_charged_particles and not contains_neutral_particles:
          self.included_particles = [3,4,5,6,7,8,9,10,11,12]
        elif contains_neutral_particles and not contains_charged_particles:
          self.included_particles = [0,1,2]
        else:
          self.included_particles = list(range(13))
        self.x_labels = particle_array[[*self.included_particles, 13]]
        self.y_labels = particle_array[self.included_particles]
      case 'all':
        self.included_particles = list(range(13))
        self.x_labels = particle_array[[*self.included_particles, 13]]
        self.y_labels = particle_array[self.included_particles]
      case 'necessary':
        self.included_particles = list(set(list(labels) + list(predictions)))
        self.included_particles.sort()
        if 13 in self.included_particles:
          self.included_particles = [i for i in self.included_particles if int(i) != 13]
        self.x_labels = particle_array[[*self.included_particles, 13]]
        self.y_labels = particle_array[self.included_particles]
      case _:
        raise ValueError("Label selection must be one of the following: 'charge','all','necessary'")

    self.nXticks = len(self.x_labels)
    self.nYticks = len(self.y_labels)

    np.add.at(temp_confusion_matrix,(labels,predictions),1)
    if np.sum(temp_confusion_matrix[:, 13]) > 0:
      temp_confusion_matrix = temp_confusion_matrix[self.included_particles, :]
      temp_confusion_matrix = temp_confusion_matrix[:,[*self.included_particles,13]]
    else:
      temp_confusion_matrix = temp_confusion_matrix[self.included_particles, :]
      temp_confusion_matrix = temp_confusion_matrix[:,self.included_particles]
      self.x_labels = [i for i in self.x_labels if i != "No ID"]
      self.nXticks -= 1

    

    if self.purity:
      temp_confusion_matrix = np.transpose(temp_confusion_matrix)
    
    self.confusion_matrix = np.zeros_like(temp_confusion_matrix)

    for i in range(len(temp_confusion_matrix)):
      self.confusion_matrix[i]= temp_confusion_matrix[i]/sum(temp_confusion_matrix[i]) if sum(temp_confusion_matrix[i]) > 0 else np.zeros_like(temp_confusion_matrix[i])

    if self.purity:
      self.confusion_matrix = np.transpose(self.confusion_matrix)

  def display_matrix(self, title):
    self.fig, self.ax = plt.subplots()
    self.im = self.ax.imshow(self.confusion_matrix)
    self.text = None
    cmap_min, cmap_max = self.im.cmap(0), self.im.cmap(1.0)
    self.text = np.empty_like(self.confusion_matrix, dtype=object)
    thresh = (self.confusion_matrix.max() + self.confusion_matrix.min()) / 2.0
    for i in range(self.nYticks):
      for j in range(self.nXticks):
        color = cmap_max if self.confusion_matrix[i][j] < thresh else cmap_min
        text_cm = round(self.confusion_matrix[i][j], 2)
        if float(text_cm) == float(0):
          text_cm = 0
        default_text_kwargs = dict(ha="center", va="center", color=color)
        text_kwargs = {**default_text_kwargs}
        self.text[i][j] = self.ax.text(j, i, text_cm, **text_kwargs)

    self.fig.colorbar(self.im, ax=self.ax)
    self.ax.set(xticks=np.arange(self.nXticks),yticks=np.arange(self.nYticks),xticklabels=[format_labels(x) for x in self.x_labels],yticklabels=[format_labels(y) for y in self.y_labels],ylabel="Generated As",xlabel="Identified As")
    self.ax.set_ylim((self.nYticks - 0.5, -0.5))
    self.fig.set_figheight(7)
    self.fig.set_figwidth(7)
    self.ax.set_title(title)
    self.fig.show()
  
  def __repr__(self):
    return "" 

def split_df(input_df, training_fraction=0.9):
  if round(training_fraction,2) == 1.:
    raise ValueError("Cannot create split dataset with such a large training fraction. Reduce the training fraction.")
  elif round(training_fraction,2) == 0:
    raise ValueError("Cannot create split dataset with such a small training fraction. Increase the training fraction.")
  if round(training_fraction,2)< 0.5:
    switch_train_test= True
    every_n_events = int(round((1-training_fraction)/training_fraction) + 1)
  else:
    switch_train_test = False
    every_n_events = int(round(training_fraction/(1-training_fraction)) + 1)
  df = input_df.copy().reset_index(drop=True)
  nRows = len(df.index)
  starting_index = 0
  training_list = []
  test_list = []
  index_list = []
  counter = 0
  while starting_index <= nRows - 1 :
    counter += 1
    ending_index = starting_index + df['Number of Hypotheses'][starting_index]
    if counter % every_n_events == 0:
      test_list.extend([not switch_train_test for _ in range(starting_index, ending_index)])
      training_list.extend([switch_train_test for _ in range(starting_index, ending_index)])
    else:
      test_list.extend([switch_train_test for _ in range(starting_index, ending_index)])
      training_list.extend([not switch_train_test for _ in range(starting_index, ending_index)])
    starting_index = ending_index
  test = df.loc[test_list].reset_index(drop=True)
  training = df.loc[training_list].reset_index(drop=True)
  return training, test

def grab_events(input_df, n_each = 5000,reverse = False, return_strings = False, allow_less=False):
  particle_list = ["Photon","KLong","Neutron","Proton","K+","Pi+","AntiMuon","Positron","AntiProton","K-","Pi-","Muon","Electron","No ID"]
  if reverse:
    df = input_df[::-1].copy().reset_index(drop=True)
  else:
    df = input_df.copy().reset_index(drop=True)
  if isinstance(df['Hypothesis'][0],str):
    df['Hypothesis']=df['Hypothesis'].apply(particle_list.index)
  if isinstance(df['Generated As'][0],str):
    df['Generated As']=df['Generated As'].apply(particle_list.index)
  nRows = len(df.index)
  starting_index = 0
  training_list = []
  include_list = []
  counter = 0
  if reverse:
    current_particle = 12
  else:
    current_particle = 0
  while starting_index <= nRows - 1 :
    if df['Generated As'][starting_index] != current_particle:
      if counter <= n_each and not allow_less:
        raise ValueError("Not enough rows in dataframe to grab " + str(n_each) + " events of " + particle_list[current_particle] + " events.")
      if reverse:
        current_particle -= 1
      else:
        current_particle += 1
      counter = 0
    counter += 1
    ending_index = starting_index + df['Number of Hypotheses'][starting_index]
    if counter <= n_each:
      include_list.extend([True for _ in range(starting_index, ending_index)])
    else:
      include_list.extend([False for _ in range(starting_index, ending_index)])
    starting_index = ending_index
  if return_strings:
    df['Hypothesis']=df['Hypothesis'].apply(lambda x: particle_list[x])
  if return_strings:
    df['Generated As']=df['Generated As'].apply(lambda x: particle_list[x])
  smaller_dataset = df.loc[include_list].reset_index(drop=True)
  return smaller_dataset

def feature_importance(model, test_data_full, target='Generated As', match_hypothesis=False,n_repetitions=3):
  test_data = grab_events(test_data_full,n_each=100, allow_less=True)
  x_test = test_data[[column for column in model.feature_names_in_]]
  y_test = test_data[target]
  importances = []
  n_features_to_shuffle = len([i for i in x_test.columns if i != "Number of Hypotheses"])
  starting_index = 0
  predictions = []
  identities = []
  hypotheses = test_data['Hypothesis'].to_list()
  length_of_df = len(x_test.index)
  dfs_to_combine = [x_test]
  total_rows = length_of_df

  for column_to_shuffle in x_test.columns:
    if column_to_shuffle == "Number of Hypotheses":
      continue
    for _ in range(n_repetitions):
      total_rows += length_of_df
      shuffled_dataframe = x_test.copy()
      shuffled_dataframe[column_to_shuffle] = shuffled_dataframe[column_to_shuffle].sample(frac=1,ignore_index=True)
      hypotheses.extend(test_data['Hypothesis'].to_list())
      dfs_to_combine.append(shuffled_dataframe)
  new_test =pd.concat(dfs_to_combine, ignore_index=True)

  index_list = []
  event_nos = []
  i = 0
  nRows = len(new_test.index)
  while starting_index < nRows - 1 :
    ending_index = starting_index + new_test['Number of Hypotheses'][starting_index]
    index_list.append((starting_index, ending_index))
    event_nos.extend([i for _ in range(starting_index, ending_index)])
    i+= 1
    starting_index = ending_index
  number_of_events = i

  if match_hypothesis:
    temp_predictions = model.predict_proba(new_test)
    new_test['Confidence'] = [max(probs) for probs in temp_predictions]
    new_test['Prediction'] = np.argmax(temp_predictions, axis=1)

    identities_grouped = new_test[['Generated As','Prediction','eventNo']].groupby('eventNo')
    identities = identities_grouped['Generated As'].head(1).to_list()

    matches_hypotheses_bool_list = (new_test['Prediction'] == new_test['Hypothesis']) | (new_test['Hypothesis'] == 13)
    matching_hypotheses = new_test.loc[matches_hypotheses_bool_list]
    grouped_df = matching_hypotheses[['Generated As','Prediction','Confidence','eventNo']].groupby('eventNo')
    max_confidence_indices = grouped_df['Confidence'].idxmax()
    predictions_temp = new_test['Prediction'].iloc[max_confidence_indices]
      
    predictions = predictions_temp.set_index('eventNo').sort_index().reindex(list(range(number_of_events)),fill_value=13).to_list
    identities = grouped_df[target].head(1).to_list()
  else:
    new_test['Prediction'] = model.predict(new_test)
    new_test['eventNo'] = event_nos
    grouped_df = new_test[['Prediction','eventNo','Generated As']].groupby('eventNo')
    predictions = grouped_df['Prediction'].agg(lambda x: x.value_counts().index[0]).to_list()
    identities = grouped_df['Generated As'].head(1).to_list()
    
  identities = np.array(identities, dtype= int)
  predictions = np.array(predictions, dtype= int)
  starting_accuracy = sklearn.metrics.accuracy_score(identities[0:1300],predictions[0:1300])
  importances = [starting_accuracy-sklearn.metrics.accuracy_score(identities[1300*(i*n_repetitions+1):1300*((i+1)*n_repetitions+1)], predictions[1300*(i*n_repetitions+1):1300*((i+1)*n_repetitions+1)]) for i in range(n_features_to_shuffle)]
  important_features = [model.feature_names_in_[i] for i in range(n_features_to_shuffle) if importances[i] > 0.005]
  important_importances = [i for i in importances if i > 0.005]
  plt.bar(important_features,important_importances)
  plt.title("Feature Importances")
  plt.xlabel("Feature Name")
  plt.ylabel("Feature Importance")
  plt.xticks(rotation=90)
  plt.show()

def save_model(model, path="my_model.joblib"):
  joblib.dump(model, path)
  print('Model saved as ' + path)

def load_model(path="my_model.joblib"):
  print("Loading model...")
  model = joblib.load(path)
  print("Done loading model")
  return model