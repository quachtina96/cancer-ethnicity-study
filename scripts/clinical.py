import numpy as np
import pickle
import sys

def get_clinical(filename):
	patients = None
	info = {}
	index = 0
	patient_info = {}

	with open(filename) as f:
		for line in f:
			# Get patient IDs
			if index == 0:
				patients = [item.upper() for item in line.split("\t")[1:]]
				end_range = len(patients) + 1

			# Get info type
			else:
				parts = line.split("\t")
				info[parts[0]] = parts[1:]
			index += 1

	assert len(info['race']) == len(patients)

	for i, patient in enumerate(patients):
		path = None

		if 'pathologic_stage' in info:
			path = info['pathologic_stage'][i]
			
		patient_info[patient] = {
			'years_to_birth': info['years_to_birth'][i],
			#'pathologic_stage': path,
			'gender': info['gender'][i],
			'race': info['race'][i],
			'ethnicity': info['ethnicity'][i]
		}
		print patient_info[patient]
	pickle.dump(patient_info, open(filename + '.saved.p','wb'))	

def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    if len(sys.argv) < 1:
        print "you must call program as: python ./clinical.py <datafile>"
        sys.exit(1)

    filename = sys.argv[1]
    get_clinical(filename)

if __name__ == "__main__":
    main()