# Imports
from ishneholterlib import Holter
import os
# len(os.listdir("./ecg_200"))
# from user_ekm_functions import user_ekm_dataset
from user_ekm_functions import user_ekm_sbf_bpf_dataset
from user_ekm_functions import user_ekm_sbf_bpf_dataset_single_thread
import multiprocessing

## Initial variables (codes for different status)

# default values for fields that may be set/edited by the user
header_field_defaults = {
    'file_version':  -9,
    'first_name':    '',
    'last_name':     '',
    'id':            '',
    'sex':           0,
    'race':          0,
    'birth_date':    None,
    'pm':            -9,
    'recorder_type': '',  # TODO?: 'unknown'
    'proprietary':   '',
    'copyright':     '',
    'reserved':      '',
    'var_block':     ''
}

# numeric codes from Table 1 of ISHNE Holter spec
lead_specs = {
    -9: 'absent', 0: 'unknown', 1: 'generic',
    2: 'X',    3: 'Y',    4: 'Z',
    5: 'I',    6: 'II',   7: 'III',
    8: 'aVR',  9: 'aVL', 10: 'aVF',
    11: 'V1', 12: 'V2',  13: 'V3',
    14: 'V4', 15: 'V5',  16: 'V6',
    17: 'ES', 18: 'AS',  19: 'AI'
}

# numeric codes from Table 2 of ISHNE Holter spec
lead_qualities = {
    -9: 'absent',
    0: 'unknown',
    1: 'good',
    2: 'intermittent noise',
    3: 'frequent noise',
    4: 'intermittent disconnect',
    5: 'frequent disconnect'
}

# type of pacemaker
pm_codes = {
    0: 'none',  # i.e. no PM installed.  so -9 should be used for unknown.
    1: 'unknown type',
    2: 'single chamber unipolar',
    3: 'dual chamber unipolar',
    4: 'single chamber bipolar',
    5: 'dual chamber bipolar',
}

gender_codes = {
    0: None,  # unknown
    1: 'M',
    2: 'F'
}

# race codes.  other values (e.g. 4+) may also be used, but weren't in the initial spec
race_codes = {
    0: None,  # unknown
    1: 'caucasian',
    2: 'black',
    3: 'oriental',
}

## Preprocessing by pan tompkins algorithm

# Load a file from disk:
all_data = Holter('./1-300m.ecg')
all_data.load_data()

x_lead = all_data.lead[0]
y_lead = all_data.lead[1]
z_lead = all_data.lead[2]

dataset_path = "./ecg_200"
users_files = os.listdir(dataset_path)
users_files.remove("clinicalData-selected")
# len(users_files)

# Getting .ecg files of users
users_ecg_files = []
for _file in users_files:
  f_extention = _file.split(".")[1]
  if f_extention == "ecg":
    users_ecg_files.append(_file)

# Initialization of dataset extracting processing
sampling_rate = all_data.sr
dataset_name = "sbf+bpf_ekm_dataset"
base_ekms_path = f'EKM_dataset'

lead_names_dict = {
    1: "x_lead",
    2: "y_lead",
    3: "z_lead"
}

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

### Multi processing

# Specify the number of processes in the pool
num_processes = multiprocessing.cpu_count()
# num_processes = 30

# Creating slices of users' ecg file for multiprocessing
# Necessety: there were some bugs in server with all amounts of users' ecg files
slices_size = num_processes
number_of_complete_slices = len(users_ecg_files)//slices_size
users_ecg_files_chunks = [users_ecg_files[_ * slices_size: (_+1) * slices_size] for _ in range(number_of_complete_slices)]

if number_of_complete_slices * slices_size != len(users_ecg_files):
    users_ecg_files_chunks.append(users_ecg_files[number_of_complete_slices * slices_size:])

def processing_ecg_files(users_ecg_files_chunk):
    with multiprocessing.Manager() as manager:
        # Create a shared counter
        shared_counter = manager.Value('i', 0)

        # Create a lock from the manager
        lock = manager.Lock()

        # Create a pool of processes
        with multiprocessing.Pool(processes=num_processes) as pool:
            # Pass the shared counter, lock, and total number of elements to the worker function
            pool.starmap(user_ekm_sbf_bpf_dataset, [(user, shared_counter, lock, len(users_ecg_files_chunk)) for user in users_ecg_files_chunk])


for users_ecg_files_chunk in users_ecg_files_chunks:
    processing_ecg_files(users_ecg_files_chunk)

# Print final progress
print("Processing complete.")

# # Single Thread EKM extractor
# for user_ecg_files in users_ecg_files:
#     user_ekm_sbf_bpf_dataset_single_thread(user_ecg_files)