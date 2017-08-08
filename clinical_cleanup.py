import sys
import os
import csv 
import optparse


# some file descriptors
ERROR_FILE = sys.stderr
OUTPUT_FILE = sys.stdout

CLIN_ATTR_COUNTS = {}

VALID_PROCESSING_TYPES = ['KEEP_ALL', 'IGNORE', 'MERGE', 'DERIVE', 'FIX_ALL', 'FIX_VALUE', 'FIX_ATTRIBUTE', 'ADD_ALL', 'ADD_GENOMIC_ALTERATIONS']
PROCESSING_TYPES_RULES = {
	'KEEP_ALL':'keep original data and clinical attribute name',
	'IGNORE':'skip attribute/data',
	'MERGE':'append data from multiple columns to normalized clinical attribute',
	'DERIVE':'derive data for normalized clinical attribute from original clinical attribute',
	'FIX_ALL':'map clinical attribute and data to normalized clinical attribute and normalized data value',
	'FIX_VALUE':'map data to normalized data value (clinical attribute stays the same)',
	'FIX_ATTRIBUTE':'map clinical attribute to normalized clinical attribute (data stays the same)',
	'ADD_ALL':'add clinical attribute and value (data did not exist in raw clinical data)',
	'ADD_GENOMIC_ALTERATIONS':'calculate the number of genomic alterations from the MAF'
}

CLINICAL_DATA_MAP = {}

CASE_ID_ATTRIBUTES = ['PATIENT_ID', 'SAMPLE_ID','OTHER_PATIENT_ID', 'OTHER_SAMPLE_ID']

NORMALIZED_ATTRIBUTE_LIST = []
POST_PROCESS_ATTRIBUTE_FILTER = []


def get_header(filename, map_clinical_data, calc_genomic_alterations):
	"""
	"""
	header = get_filtered_header(get_file_header(filename))
	if map_clinical_data:
		header = get_normalized_header(header)
	if calc_genomic_alterations:
		header.append('GENOMIC_ALTERATIONS')
		print 'Added GENOMIC_ALTERATIONS to clinical attributes for header.'
	return header


def get_file_header(filename):
	"""
	"""	
	filedata = open(filename).read().split('\n')	
	file_header = map(str.strip, filedata[0].split('\t'))

	header = [hdr for hdr in CASE_ID_ATTRIBUTES if hdr in file_header]
	header_ext = [hdr for hdr in file_header if hdr not in header]

	header.extend(header_ext)
	return header


def get_filtered_header(header):
	"""
	"""	
	new_header = [hdr for hdr in header if CLIN_ATTR_COUNTS.get(hdr,0) > 0]
	difference = len(header)-len(new_header)
	if difference == 0:
		print 'No clinical attributes were found with zero counts.'
	else:
		print 'Filtering out', str(difference)+'/'+str(len(header)), 'clinical attributes with zero counts.'
	return new_header


def get_normalized_header(header):
	"""
	"""
	new_header = [hdr for hdr in header if hdr in NORMALIZED_ATTRIBUTE_LIST]
	header_ext = [attr for attr in NORMALIZED_ATTRIBUTE_LIST if attr not in new_header]

	print 'Removed', str(len(header) - len(new_header)), 'clinical attributes from filtered header.'
	print 'Added', str(len(header_ext)), 'normalized clinical attributes to final header.'

	new_header.extend(header_ext)
	return new_header


def process_datum(val):
	"""
	"""
	try: 
		vfixed = val.strip()
	except AttributeError:
		vfixed = 'NA'

	if vfixed in ['', None, 'N/A']:
		return 'NA'
	else:
		return vfixed


def update_attribute_counts(line_data):
	"""
	"""
	for k,v in line_data.items():
		if v != 'NA':
			CLIN_ATTR_COUNTS[k] = CLIN_ATTR_COUNTS.get(k, 0) + 1


def basic_clinical_cleanup(clin_filename):
	"""
	"""
	clin_file = open(clin_filename, 'rU')
	clin_reader = csv.DictReader(clin_file, dialect = 'excel-tab')
	clin_attrs = get_file_header(clin_filename)

	basic_sample_data = {}
	num_vals_changed = 0
	oncotre_counts = 0
	for line in clin_reader:
		sample_id = line['SAMPLE_ID'].strip()
		line_data = map(lambda x: process_datum(line.get(x,'NA')), clin_attrs)
		for i,attr in enumerate(clin_attrs):
			if line_data[i] != line.get(attr,'NA'):
				num_vals_changed += 1
				print 'Value fixed in column:', attr, 'for sample id', sample_id
				print '\t', '"'+line.get(attr,'NA')+'"', '==>', '"'+line_data[i]+'"'
		sample_data = dict(zip(clin_attrs, line_data))
		basic_sample_data[sample_id] = sample_data
		update_attribute_counts(sample_data)

	if num_vals_changed == 0:
		print 'No values were fixed.'
	else:
		print 'Fixed,', num_vals_changed, 'values.'

	if oncotre_counts > 0:
		print 'ONCOTREE_CODE added to', oncotre_counts, 'samples'

	clin_file.close()
	return basic_sample_data


def write_temp_file(filename):
	"""
	"""
	filedata = [x for x in open(filename).read().split('\n') if not x.startswith('#')]
	temp_filename = os.path.join(os.path.dirname(filename), 'temp_'+os.path.basename(filename))
	fh = open(temp_filename, 'w')
	fh.write('\n'.join(filedata))
	fh.close()
	return temp_filename
	

def calculate_genomic_alterations(maf_filename):
	"""
	"""
	temp_filename = write_temp_file(maf_filename)
	print 'Writing temp MAF file to:', temp_filename
	
	maf_file = open(temp_filename, 'rU')
	maf_reader = csv.DictReader(maf_file, dialect = 'excel-tab')

	genomic_alts_data = {}
	for line in maf_reader:
		sample_id = line['Tumor_Sample_Barcode'].strip()
		genomic_alts_data[sample_id] = genomic_alts_data.get(sample_id, 0) + 1

	maf_file.close()	
	os.remove(temp_filename)
	print 'Removed temp MAF file.'	
	return genomic_alts_data


def normalize_attribute_data(ptype, norm_attr, sample_data):
	"""
	"""
	norm_attr_data = CLINICAL_DATA_MAP[ptype].get(norm_attr)
	orig_attrs = list(norm_attr_data.keys())

	orig_vals = map(lambda x: sample_data.get(x, 'NA'), orig_attrs)
	
	if ptype == 'MERGE':
		value = []
		for val in orig_vals:
			value.extend(val.split('/'))

		if not value:
			norm_val = 'NA'
		else:
			norm_val_list = []
			for i,val in enumerate(value):
				new_val = norm_attr_data[orig_attrs[i]].get(val, 'NA')
				if new_val == None:
					continue
				norm_val_list.extend(new_val.split('/'))

			norm_val_list = list(set(norm_val_list))
			if len(norm_val_list) > 0:
				norm_val_list = [v for v in norm_val_list if v != 'None']
				
			if not norm_val_list:
				norm_val = 'NA'
			else:
				norm_val = '/'.join(sorted(norm_val_list))

	else:
		if len(set(orig_vals)) > 1:
			print 'ERROR: more than one original value retrieved from matched original attribute.'
			print 'PTYPE:', ptype
			print 'NORM_ATTRIBUTE:', norm_attr
			print 'NORM_ATTRIBUTE_DATA:', norm_attr_data
			print 'ORIGINAL_ATTRIBUTE(S):', orig_attrs
			print 'ORIGINAL_VALUE(S):', orig_vals
			sys.exit(2)
		else:
			if orig_vals[0] == 'NA':
				norm_val = orig_vals[0]
			else:
				try:
					orig_attr_data = norm_attr_data.get(orig_attrs[0])
					norm_val = orig_attr_data.get(orig_vals[0])
					if norm_val == None:
						if orig_vals[0] in orig_attr_data.values():
							norm_val = orig_vals[0]
					
				except AttributeError:
					print 'EXCEPTION: AttributeError'
					print 'NORM_ATTR_DATA:', norm_attr_data
					print 'ORIG_ATTRS:', orig_attrs
					print 'ORIG_VALS:', orig_vals
					print '\n\n'					
	return norm_val


def get_processing_type(norm_attr):
	"""
	"""
	for ptype in CLINICAL_DATA_MAP.keys():
		ptype_data = CLINICAL_DATA_MAP.get(ptype)

		if ptype in ['KEEP_ALL', 'IGNORE']:
			if norm_attr in ptype_data:
				return ptype		
		elif ptype in ['MERGE', 'DERIVE', 'FIX_ALL', 'FIX_VALUE', 'FIX_ATTRIBUTE', 'ADD_ALL']:
			if norm_attr in ptype_data.keys():
				return ptype

	print 'ERROR: Could not find processing type for normalized attribute:', norm_attr
	return None


def get_normalized_sample_data(sample_data, header):
	"""
	"""
	for norm_attr in header:
		if norm_attr in CLINICAL_DATA_MAP.get('KEEP_ALL'):
			continue
		elif norm_attr == 'GENOMIC_ALTERATIONS':
			continue

		ptype = get_processing_type(norm_attr)
		if ptype in ['MERGE', 'DERIVE', 'FIX_ALL', 'FIX_VALUE', 'FIX_ATTRIBUTE']:
			norm_val = normalize_attribute_data(ptype, norm_attr, sample_data)
			sample_data[norm_attr] = norm_val
		elif ptype == 'ADD_ALL':
			norm_val = CLINICAL_DATA_MAP[ptype].get(norm_attr)
			sample_data[norm_attr] = norm_val
		else:
			print 'ERROR: Attribute in header has not been normalized.'
			print norm_attr
			sys.exit(2)

	processed_sample_data = map(lambda x: sample_data.get(x, 'NA'), header)
	return processed_sample_data


def cleanup_clinical_data(clin_filename, output_directory, map_clinical_data, calc_genomic_alterations):
	"""
	"""
	basic_sample_data = basic_clinical_cleanup(clin_filename)
	header = get_header(clin_filename, map_clinical_data, calc_genomic_alterations)
	# if 'ONCOTREE_CODE' not in header:
	# 	header.append('ONCOTREE_CODE')

	if calc_genomic_alterations:
		maf_filename = os.path.join(os.path.dirname(clin_filename), 'data_mutations_extended.txt')
		genomic_alts_data = calculate_genomic_alterations(maf_filename)

	output = ['\t'.join(header)]
	for sample_id,sample_data in basic_sample_data.items():
		if calc_genomic_alterations:
			sample_data['GENOMIC_ALTERATIONS'] = str(genomic_alts_data.get(sample_id, 0))

		if map_clinical_data:
			processed_sample_data = get_normalized_sample_data(sample_data, header)
		else:
			processed_sample_data = map(lambda x: sample_data.get(x, 'NA'), header)

		output.append('\t'.join(processed_sample_data))

	output_filename = os.path.join(output_directory, 'processed-'+os.path.basename(clin_filename))
	fh = open(output_filename, 'w')
	fh.write('\n'.join(output))
	fh.close()

	print 'Filtered clinical data written to:', os.path.abspath(output_filename)


def generate_clinical_data_map(filename):
	"""
	"""
	map_file = open(filename, 'rU')
	map_reader = csv.DictReader(map_file, dialect = 'excel-tab')

	for line in map_reader:
		ptype = line.get('PROCESSING_TYPE','NA').strip()
		if ptype not in VALID_PROCESSING_TYPES:
			print 'Invalid processing type:', ptype
			processing_type_rules()
			sys.exit(2)

		orig_attr = line.get('ORIGINAL_ATTRIBUTE', 'NA').strip()
		orig_val = line.get('ORIGINAL_VALUE', 'NA').strip()		
		norm_attr = line.get('NORMALIZED_ATTRIBUTE', 'NA').strip()
		norm_val = line.get('NORMALIZED_VALUE', 'NA').strip()						

		if ptype in ['KEEP_ALL', 'IGNORE']:
			try:
				CLINICAL_DATA_MAP[ptype].append(orig_attr)
			except KeyError:
				CLINICAL_DATA_MAP[ptype] = [orig_attr]

		elif ptype == 'ADD_ALL':
			ptype_data = CLINICAL_DATA_MAP.get(ptype, {})
			ptype_data[norm_attr] = norm_val
			CLINICAL_DATA_MAP[ptype] = ptype_data

		elif ptype in ['MERGE', 'DERIVE', 'FIX_ALL', 'FIX_VALUE', 'FIX_ATTRIBUTE']:
			ptype_data = CLINICAL_DATA_MAP.get(ptype, {})
			norm_attr_data = ptype_data.get(norm_attr, {})
			orig_attr_data = norm_attr_data.get(orig_attr, {})

			orig_attr_data[orig_val] = norm_val
			norm_attr_data[orig_attr] = orig_attr_data
			ptype_data[norm_attr] = norm_attr_data
			CLINICAL_DATA_MAP[ptype] = ptype_data		

		# add normalized clinical attributes to NORMALIZED_ATTRIBUTE_LIST if not already in list
		if ptype in ['KEEP_ALL', 'ADD_ALL', 'DERIVE', 'FIX_VALUE', 'FIX_ATTRIBUTE' ,'FIX_ALL', 'MERGE']:
			if norm_attr not in NORMALIZED_ATTRIBUTE_LIST:
				NORMALIZED_ATTRIBUTE_LIST.append(norm_attr)
		
		# add original clinical attributes to POST_PROCESS_ATTRIBUTE_FILTER if meets conditions below
		if ptype in ['MERGE', 'FIX_ALL'] and orig_attr != norm_attr and orig_attr not in POST_PROCESS_ATTRIBUTE_FILTER:
				POST_PROCESS_ATTRIBUTE_FILTER.append(orig_attr)
		elif ptype == 'IGNORE' and orig_attr not in POST_PROCESS_ATTRIBUTE_FILTER:
			POST_PROCESS_ATTRIBUTE_FILTER.append(orig_attr)
	print CLINICAL_DATA_MAP['MERGE']


def processing_type_rules():
	print
	print '\tPROCESSING_TYPE\tRULE'
	print '\t----------------------'
	for ptype in PROCESSING_TYPES:
		print '\t'+ptype+'\t'+PROCESSING_TYPES_RULES.get(ptype)


def usage():
	print >> OUTPUT_FILE, 'clinical_cleanup.py --clinical-file [path/to/clinical/file] --output-directory [path/to/output/directory] --map-file [path/to/map/file] --genomic-alterations [True/False]'


def main():
	# get command line stuff
	parser = optparse.OptionParser()
	parser.add_option('-c', '--clinical-file', action = 'store', dest = 'clinfile')
	parser.add_option('-d', '--output-directory', action = 'store', dest = 'outputdir')
	parser.add_option('-m', '--map-file', action = 'store', dest = 'mapfile')
	parser.add_option('-g', '--genomic-alterations', action = 'store', dest = 'genalts')

	(options, args) = parser.parse_args()

	clin_filename = options.clinfile
	output_directory = options.outputdir
	map_filename = options.mapfile
	genomic_alterations = options.genalts


	map_clinical_data = False
	calc_genomic_alterations = False


	if clin_filename == None:
		usage()
		sys.exit(2)
	elif not os.path.exists(clin_filename):
		print 'ERROR: No such file:', os.path.abspath(clin_filename)
		sys.exit(2)

	if output_directory == None:
		print 'WARNING: No output directory defined - using clinical file directory instead.'
		output_directory = os.path.dirname(clin_filename)

	if map_filename == None:
		print 'WARNING: No clinical data map provided.'
		print 'WARNING: Will not normalize clinical data.'
	else: 
		if os.path.exists(map_filename):			
			print 'Mapping clinical attributes and data using:', os.path.abspath(map_filename)
			generate_clinical_data_map(map_filename)
			map_clinical_data = True
		else:
			print 'ERROR: No such file:', map_filename
			sys.exit(2)

	
	if genomic_alterations == None:
		print 'WARNING: No genomic alterations argument found - will not add to clinical file.'
	else:
		genomic_alterations = genomic_alterations.strip().lower()
		if genomic_alterations not in ['true', 'false']:
			print 'ERROR: Invalid argument for genomic alterations:', genomic_alterations
			usage()
			sys.exit(2)				
		elif genomic_alterations == 'false':
			print 'Genomic alterations will not be calculated for clinical file.'
		elif genomic_alterations == 'true':
			maf_filename = os.path.join(os.path.dirname(clin_filename), 'data_mutations_extended.txt')
			if not os.path.exists(maf_filename):
				print 'ERROR: MAF file could not be found in directory:', os.path.dirname(clin_filename)
				print 'ERROR: Please make sure that MAF file exists in same directory as clinical file.'
				sys.exit(2)
			else:
				print 'Genomic alterations will be calculated from:', os.path.abspath(maf_filename)
				calc_genomic_alterations = True

	cleanup_clinical_data(clin_filename, output_directory, map_clinical_data, calc_genomic_alterations)


if __name__ == '__main__':
	main()