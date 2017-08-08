import sys
import os
import argparse
import re
import xml.etree.ElementTree as ET
from StringIO import StringIO


def fmi_xml_processor(filename,output_directory):

	print 'Processing filename:', filename

	tree = ET.parse(filename)
	root = tree.getroot()
	variant_report_list = root.getchildren()
	
	new_tree = generate_multi_case_tree(variant_report_list)

	output_filename = os.path.join(output_directory, os.path.basename(filename))
	output_file = ET.ElementTree(new_tree)
	output_file.write(output_filename)

	print 'Processed xml written to:', output_filename

def generate_multi_case_tree(variant_report_list):
	tree = ET.Element('ClientCaseInfo')
	root = ET.SubElement(tree, 'Cases')

	for variant_report in variant_report_list:
		case,data = get_case_data(variant_report)
		variant_report.attrib['test-request'] = case
		variant_report = process_variant_report(root, variant_report)
		case_element = ET.SubElement(root, 'Case', data)
		case_element.append(variant_report)

	return tree

def process_variant_report(root, variant_report):
	if not 'non-human-content' in variant_report.getchildren():
		print 'Adding non-human content'
		non_human_content = ET.SubElement(root, 'non-human-content', {})
		variant_report.append(non_human_content)
	if not 'pipeline-version' in variant_report.keys():
		print 'Adding pipeline version'
		variant_report.attrib['pipeline-version'] = ''

	short_variants = variant_report.find('short-variants')
	for i,short_variant in enumerate(variant_report.find('short-variants').getchildren()):
		new_short_variant = short_variant
		if not 'depth' in short_variant.keys():
			print 'Adding depth to short variant'
			new_short_variant.attrib['depth'] = ''
		variant_report.find('short-variants').remove(short_variant)
		variant_report.find('short-variants').append(new_short_variant)
	return variant_report


def get_case_data(variant_report):
	try:
		case = variant_report.find('samples').find('sample').get('name')
	except:
		print 'Using specimen - could not find samples element in variant-report'
		case = variant_report.get('specimen')
	data = {'case':case, 'fmiCase':case, 'hasVariant':'1'}
	return case,data


def get_data_filenames(data_directory):
	filenames = [f for f in os.listdir(data_directory) if f.endswith('.xml')]

	output_directory = os.path.join(os.path.dirname(datasource),'stripped')
	if not os.path.exists(output_directory):
		print 'Creating output directory:', output_directory
		os.makedirs(output_directory)

	case_filenames = {}
	datadir = os.path.join(os.path.dirname(datasource), 'orig')
	for fname in filenames:
		case = os.path.basename(fname).split('.')[0]
		case_filenames[case] = os.path.join(datadir, os.path.basename(fname))

	return case_filenames,output_directory


def cleanup_xml_file(output_directory, filename):
	output_filename = os.path.join(output_directory, filename)
	tree = ET.parse(filename)
	root = tree.getroot()

	it = ET.iterparse(StringIO(ET.tostring(root)))
	for _, el in it:
		if '}' in el.tag:
			el.tag = el.tag.split('}', 1)[1]  # strip all namespaces
	root = it.root	

	ntree = ET.Element('ClientCaseInfo')
	ntree.extend(root.getchildren())
	output = ET.ElementTree(ntree)
	output.write(output_filename)
	print 'File written to: ' + output_filename
	return output_filename


def interface(args=None):
	parser = argparse.ArgumentParser(description='Foundation XML data cleanup.'
    								'\noptions must include either --data-directory or --filename, but not both')
	parser.add_argument('-d', '--data_directory', type=str, required=False, help='Path to XML data directory')
	parser.add_argument('-f', '--filename', type=str, required=False, help='Path to XML file')
	parser.add_argument('-o', '--output_directory', type=str, required=True, help='Path to output directory for stripped XML files')
	return parser.parse_args(args)


def main():
	parsed_args = interface()

	data_directory = parsed_args.data_directory
	output_directory = parsed_args.output_directory
	filename = parsed_args.filename

	files = []
	if filename != None and os.path.exists(filename):
		new_filename = cleanup_xml_file(output_directory, filename)
		fmi_xml_processor(new_filename, output_directory)
	elif data_directory != None and os.path.exists(data_directory):
		for fname in os.path.listdir(data_directory):
			filename = os.path.join(data_directory, fname)
			new_filename = cleanup_xml_file(output_directory, filename)
			fmi_xml_processor(new_filename, output_directory)

if __name__ == '__main__':
	main()