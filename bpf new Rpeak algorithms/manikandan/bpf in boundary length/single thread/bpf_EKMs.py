########################################
#       Installation and Imports
########################################

import os
from user_ekms_functions import user_ekm_dataset

########################################
#           Initial variables
########################################

dataset_path = "../../../../datasets/ECG 200 dataset/ecg200"
users_files = os.listdir(dataset_path)
users_files.remove("clinicalData-selected")

# Initialization of dataset extracting processing
dataset_name = "main_ekm_dataset"
base_ekms_path = f'EKM_dataset'

lead_names_dict = {
    1: "x_lead",
    2: "y_lead",
    3: "z_lead"
}

if not os.path.isdir('./Users EKM zip'):
  os.mkdir("./Users EKM zip")

########################################
#       Reading files and playground!
########################################

# Getting .ecg files of users
users_ecg_files = []
for _file in users_files:
  f_extention = _file.split(".")[1]
  if f_extention == "ecg":
    users_ecg_files.append(_file)

# Removing users that their EKMs have been extracted
user_ekms_zip_files = os.listdir("./Users EKM zip")
list_of_ekm_extracted_users = []

for _file in user_ekms_zip_files:
  try:
    extention = _file.split(".")[1]
    if extention == "zip":
      user_id = _file.split(".")[0]
      list_of_ekm_extracted_users.append(user_id)
  except:
    pass

for usr in list_of_ekm_extracted_users:
  try:
    users_ecg_files.remove(usr + ".ecg")
  except:
    print(f"User No. {usr} already been removed.")

# Processing the ecg files
for _, ecg_file in enumerate(users_ecg_files):
    print(ecg_file)
    user_ekm_dataset(ecg_file)
    
    if _ % 10 == 0:
      print(f"Progress: {(_+1)/len(users_ecg_files)}%")

# Print final progress
print("Processing complete.")