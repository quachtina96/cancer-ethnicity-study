'''
@fileoverview This file provides functions for getting demographic information, 
which is necessary for our cancer ethnicity study but is not associated with the 
samples in which SNPs were found in the Single Nucleotide Variation MAF files. 

@author Tina Quach (quacht@mit.edu)
'''

import os
from subprocess import call
import json

# Given a tumor_sample_barcode string, queries the Genomic Data Commons through their
# API in order to retrieve the Case Id associated with the tumor sample barcode.
def getCaseId(tumor_sample_barcode):
	file = tumor_sample_barcode + '_case_id.txt'

	# Form the command (that you would normally type into the command line) to
	# query the database and write the JSON response to a file
	get_case_id = "curl -s 'https://gdc-api.nci.nih.gov/files/ids?query=" + tumor_sample_barcode + "&pretty=true' -o " + file
	os.system(get_case_id)

	# Parse the JSON response
	case_id_query_output = open(file, 'r').read()
	try:
		case_json = json.loads(case_id_query_output)
	except: 
		print "Unexpected error: Could not load JSON from file: " + file
	
	try:
		case_id = case_json["data"]["hits"][0]["cases"][0]["case_id"]
	except:
		print "Unexpected error: Could not find case ID corresponding to Tumor Sample Barcode: %s" %(tumor_sample_barcode)
		case_id = None
	os.remove(file)
	return case_id

def createCasePayload(caseId):
	# Create a payload
	payload_name = caseId+'_payload'
	payload = open(payload_name, 'w')
	payload_content = '''{
    "filters":{
        "op":"in",
        "content":{
            "field":"cases.case_id",
            "value":[
                "''' + caseId + '''"
            ]
        }
    },
    "format":"JSON",
    "fields":"case_id,demographic.demographic_id,demographic.ethnicity,demographic.race",
    "size":"100"
	}'''
	payload.write(payload_content)
	return payload_name

def getDemographicInfo(caseId):
	payload_name = createCasePayload(caseId)
	cases_endpoint = '"https://gdc-api.nci.nih.gov/cases"'
	file = caseId +'_demographic.txt'

	get_demographic = 'curl -s --request POST --header "Content-Type: application/json" --data @' + payload_name + ' ' + cases_endpoint + ' | python -m json.tool > ' + file 
	os.system(get_demographic)

	# os.system('cat temp_demographic.txt')
	# Parse the JSON response
	demographic_query_output = open(file, 'r').read()
	try:
		demographic_json = json.loads(demographic_query_output)
	except: 
		print "Unexpected error: Could not load JSON from file: " + file
	
	try:
		demographic_dict = demographic_json["data"]["hits"][0]["demographic"]
		del demographic_dict["demographic_id"]
	except:
		print "Unexpected error: Could not find demographic information for case."
	os.remove(file)
	return demographic_dict

def createTumorUUIDPayload(tumorUUID):
	# Create a payload
	payload_name = tumorUUID+'_payload'
	payload = open(payload_name, 'w')
	payload_content = '''{
    "filters":{
        "op":"in",
        "content":{
            "field":"cases.aliquot_ids",
            "value":[
                "''' + tumorUUID + '''"
            ]
        }
    },
    "format":"JSON",
    "fields":"case_id,demographic.demographic_id,demographic.ethnicity,demographic.race",
    "size":"100"
	}'''
	payload.write(payload_content)
	return payload_name

def getDemographicFromTumorUuid(tumorUUID):
	payload_name = createTumorUUIDPayload(tumorUUID)
	cases_endpoint = '"https://gdc-api.nci.nih.gov/cases"'
	file = tumorUUID + '_demographic.txt'

	get_demographic = 'curl -s --request POST --header "Content-Type: application/json" --data @' + payload_name + ' ' + cases_endpoint + ' | python -m json.tool > ' + file 
	os.system(get_demographic)

	# Parse the JSON response
	demographic_query_output = open(file, 'r').read()
	demographic_dict = None
	try:
		demographic_json = json.loads(demographic_query_output)
	except: 
		print "Unexpected error: Could not load JSON from file: " + file
	
	try:
		demographic_dict = demographic_json["data"]["hits"][0]["demographic"]
		del demographic_dict["demographic_id"]
	except:
		print "Unexpected error: Could not find demographic information for the case associated with the Tumor Sample UUID: %s" %(tumorUUID)
	os.remove(file)
	os.remove(payload_name)
	return demographic_dict


if __name__ == '__main__':
	tumor_sample_barcode = 'TCGA-3C-AALJ-01A-31D-A41F-37925'
	tumorUUID = 'c13e05b6-0b3e-4c1f-800e-4a49c5d417aa'
	assert({u'race': u'white', u'ethnicity': u'not hispanic or latino'} == getDemographicFromTumorUuid(tumorUUID))
	



