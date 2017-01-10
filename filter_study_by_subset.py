import os
import sys
import optparse
import csv

PATIENT_ID_LIST = []
SAMPLE_SUBSET_LIST = []

NON_CASE_IDS = ['Hugo_Symbol', 'Entrez_Gene_Id']
PROFILE_DATATYPE_FILENAMES = [
	'data_CNA.txt', 'data_expression_median.txt', 'data_expression_miRNA.txt', 
	'data_methylation_hm27.txt', 'data_RNA_Seq_expression_median.txt'
]

NORMAL_DATATYPE_FILENAMES = {
	'data_clinical.txt':'SAMPLE_ID',
	'data_clinical_patient.txt':'PATIENT_ID',
	'data_clinical_sample.txt':'SAMPLE_ID',
	'data_mutations_extended.txt':'Tumor_Sample_Barcode',
	'_data_cna_hg18.seg':'ID',
	'_data_cna_hg19.seg':'ID',
	'data_timeline':'PATIENT_ID_LIST'
}


def process_datum(value):
	""" 
		Strips leading/trailing whitespace from datum. 
		Returns empty string if AttributeError exception is caught. 
	"""
	try: 
		vfixed = value.strip()
	except AttributeError:
		vfixed = ''
	return vfixed


def get_header(filename):
	""" Gets the header from the file. """

	filedata = [x for x in open(filename).read().split('\n') if not x.startswith("#")]
	header = map(str.strip, filedata[0].split('\t'))
	return header


def filter_profile_data_file(filename, subset_id):
	""" Filters profile data file sample subset list """

	# flag for indicating whether data was successfully filtered or not
	filter_sucessful = True

	# get header for subset data
	subset_header = [hdr for hdr in get_header(filename) if hdr in NON_CASE_IDS]
	non_case_id_count = len(subset_header)
	subset_header.extend(SAMPLE_SUBSET_LIST[:])

	# open data file and load data only for samples in subset list
	data_file = open(filename, 'rU')
	data_reader = csv.DictReader(data_file, dialect='excel-tab')

	filtered_data = ['\t'.join(subset_header)]
	for i,line in enumerate(data_reader):
		# read data only for columns in new subset header
		row_data = map(lambda x: line.get(x,''), subset_header)

		# count number of empty columns
		non_empty_data = [val for val in row_data if val != '']

		# make sure that row data is not empty, exit if true
		if len(non_empty_data) <= non_case_id_count:
			filter_sucessful = False
			break

		# clean up data values before concatenating and appending to filtered data list
		processed_row_data = map(lambda x: process_datum(x), row_data)
		filtered_data.append('\t'.join(processed_row_data))
	data_file.close()

	# if filtering not successful then alert
	if not filter_sucessful:
		print 'ERROR: Data could not be filtered using given subset list - skipping file:', filename
		print
	else:
		# write filtered data to file 
		print 'Data successfully filtered for file:', filename
		output_filename = filename + '.' + subset_id
		output_file = open(output_filename, 'w')
		output_file.write('\n'.join(filtered_data))
		output_file.close()

		print 'Filtered data written to:', output_filename
		print


def filter_normal_data_file(filename, column, subset_id):
	""" Filters data file by sample subset list using given case id column. """

	# get header for subset data
	header = get_header(filename)

	# open data file and load data only for samples in subset list
	data_file = open(filename, 'rU')
	data_reader = csv.DictReader(data_file, dialect='excel-tab')

	filtered_data = ['\t'.join(header)]
	for i,line in enumerate(data_reader):
		# skip row of data if sample id not in sample subset list
		if line[column] in SAMPLE_SUBSET_LIST or line[column] in PATIENT_ID_LIST:			
			row_data = map(lambda x: line.get(x,''), header)

			# clean up data values before concatenating and appending to filtered data list
			processed_row_data = map(lambda x: process_datum(x), row_data)
			filtered_data.append('\t'.join(processed_row_data))
	data_file.close()

	# if filtering not successful then alert
	if len(filtered_data) == 1:
		print 'ERROR: Data could not be filtered using given subset list - skipping file:', filename
		print
	else:
		# write filtered data to file 
		print 'Data successfully filtered for file:', filename
		output_filename = filename + '.' + subset_id
		output_file = open(output_filename, 'w')
		output_file.write('\n'.join(filtered_data))
		output_file.close()

		print 'Filtered data written to:', output_filename
		print


def filter_study_by_subset_main(input_directory, subset_id):
	""" Filters study by sample subset list provided. """
	for filename in os.listdir(input_directory):
		# skip sub-directories
		if os.path.isdir(filename) or 'meta_' in filename:
			print 'Skipping sub-directory or meta file:', filename
			print
			continue

		if len(SAMPLE_SUBSET_LIST) == 0:
			print 'ERROR: sample subset list got deleted'
			sys.exit(2)

		filename_path = os.path.join(input_directory, filename)
		if filename in PROFILE_DATATYPE_FILENAMES:
			print 'Processing data from file:', filename
			filter_profile_data_file(filename_path, subset_id)

		elif filename in NORMAL_DATATYPE_FILENAMES.keys() or filename.endswith('.seg'):
			# get case id column name
			if filename.endswith('.seg'):
				column = 'ID'
			else:
				column = NORMAL_DATATYPE_FILENAMES[filename]

			print 'Processing data from file:', filename, 'using column:', column
			filter_normal_data_file(filename_path, column, subset_id)
		else:
			print 'Skipping unknown filename pattern:', filename
		print


def generate_patient_subset_list(input_directory):
	""" Generates list of patients using  """

	# determine which clinical filename to use for generating list of patients from sample subset list
	clinical_filenames = [filename for filename in os.listdir(input_directory) if 'data_clinical' in filename]
	if len(clinical_filenames) == 0:
		print 'ERROR: No clinical data files found in:', input_directory
		sys.exit(2)
	if 'data_clinical.txt' in clinical_filenames:
		clinical_filename = os.path.join(input_directory, 'data_clinical.txt')
	else:
		clinical_filename = os.path.join(input_directory, 'data_clinical_sample.txt')

	# open clinical file and generate list of patient ids
	clin_file = open(clinical_filename, 'rU')
	clin_reader = csv.DictReader(clin_file, dialect='excel-tab')
	for line in clin_reader:
		if line['SAMPLE_ID'] in SAMPLE_SUBSET_LIST and not line['PATIENT_ID'] in PATIENT_ID_LIST:
			PATIENT_ID_LIST.append(line['PATIENT_ID'])
	clin_file.close()


def load_sample_subset_list(input_directory, subset_filename):
	""" Loads subset list from given file. """

	subset_file = open(subset_filename, 'rU')	
	SAMPLE_SUBSET_LIST.extend(map(str.strip, subset_file.read().split('\n')))
	subset_file.close()

	# generate patient subset list 
	generate_patient_subset_list(input_directory)

def usage():
	print >> OUTPUT_FILE, 'filter_study_by_subset.py --subset-file path/to/subset/file --subset-identifier identifier_for_subset [default=filtered] --input-directory path/to/input/directory'

def main():
	parser = optparse.OptionParser()
	parser.add_option('-f', '--subset-file', action = 'store', dest = 'subsetfile')
	parser.add_option('-s', '--subset-identifier', action='store', dest='subsetid')
	parser.add_option('-i', '--input-directory', action='store', dest='inputdir')

	(options, args) = parser.parse_args()

	subset_filename = options.subsetfile
	subset_id = options.subsetid
	input_directory = options.inputdir

	if not os.path.exists(subset_filename):
		print 'No such file or directory:', subset_filename
		sys.exit(2)

	if not os.path.exists(input_directory):
		print 'No such file or directory:', input_directory
		sys.exit(2)

	if not subset_id:
		print 'No subset identifier entered - using default value "filtered"'

	# load sample subset list from file
	load_sample_subset_list(input_directory, subset_filename)

	# filter study by samples in subset list
	filter_study_by_subset_main(input_directory, subset_id)


if __name__ == '__main__':
	main()