import os
import sys
import csv
import optparse

# some file descriptors
ERROR_FILE = sys.stderr
OUTPUT_FILE = sys.stdout

# file types
MUTATION_FILE_PATTERN = 'data_mutations_extended.txt'
CLINICAL_FILE_PATTERN = 'data_clinical.txt'
CLINICAL_PATIENT_FILE_PATTERN = 'data_clinical_patient.txt'
CLINICAL_SAMPLE_FILE_PATTERN = 'data_clinical_sample.txt'

CASE_ID_MAP = {MUTATION_FILE_PATTERN:'Tumor_Sample_Barcode',
	CLINICAL_FILE_PATTERN:'SAMPLE_ID',
	CLINICAL_PATIENT_FILE_PATTERN:'SAMPLE_ID',
	CLINICAL_SAMPLE_FILE_PATTERN:'SAMPLE_ID'
}

def get_header(filename):
	""" Returns the file header. """	
	filedata = [x for x in open(filename).read().split('\n') if not x.startswith('#')]
	header = map(str.strip, filedata[0].split('\t'))
	return header

def get_case_ids(filename, column):
	""" Returns the case ids from the file. """
	cases = []

	# get the case ids from the designated column in the file
	data_file = open(filename, 'rU')
	file_reader = csv.DictReader(data_file, dialect='excel-tab')
	for line in file_reader:
		case_id = line.get(column, 'NA').strip()
		cases.append(case_id)
	data_file.close()

	return list(set(cases))


def insert_sequenced_samples_main(source_file, maf_file, output_directory):
	""" Creates a new MAF with the sequenced samples tag. """

	# get the case id column from the source file and create list of case ids 
	id_column = CASE_ID_MAP.get(os.path.basename(source_file), 'SAMPLE_ID')
	case_id_list = get_case_ids(source_file, id_column)

	# get everything from the file except any commented lines	
	filedata = [x for x in open(maf_file).read().split('\n') if not x.startswith('#')]
	sequenced_samples_tag = '#sequenced_samples: ' + ' '.join(case_id_list)

	# format the data for the output file
	output_data = '\n'.join([sequenced_samples_tag, '\n'.join(filedata)])

	output_filename = os.path.join(output_directory, 'data_mutations_extended_seqsamples.txt')
	output_file = open(output_filename, 'w')
	output_file.write(output_data)
	output_file.close()
	print 'MAF with sequenced samples tag written to:', output_filename


def usage():
	print >> OUTPUT_FILE, 'insert_sequenced_samples.py --source-file path/to/source --maf-file [path/to/maf] --output-directory [/path/to/output]'


def main():
	# get command line arguments
	parser = optparse.OptionParser()
	parser.add_option('-s', '--source-file', action = 'store', dest = 'sourcefile')
	parser.add_option('-d', '--output-directory', action = 'store', dest = 'outputdir')
	parser.add_option('-m', '--maf-file', action = 'store', dest = 'maffile')

	(options, args) = parser.parse_args()
	source_file = options.sourcefile
	# output_directory = options.outputdir
	output_directory = os.path.dirname(source_file)
	maf_filename = options.maffile

	# set output directory to source file directory if none provided
	if not output_directory:
		output_directory = os.getcwd()
		

	# exit it source file does not exist
	if not os.path.exists(source_file):
		print 'No such file:', source_file
		sys.exit(2)

	# exit if output directory does not exist
	if not os.path.exists(output_directory):
		print 'No such directory:', output_directory
		sys.exit(2)

	insert_sequenced_samples_main(source_file, maf_filename, output_directory)


if __name__ == '__main__':
	main()
