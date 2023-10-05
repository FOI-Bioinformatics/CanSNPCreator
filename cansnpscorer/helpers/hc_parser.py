import csv
import sys, re

def clean_name(name):
	return re.sub('[^A-Za-z0-9]+', '', name)

class hc_parser():
	def __init__(self,csvfile):
		self.dictonary = {}
		self.group_label = "Group"
		self.group_index = 0
		self.strain_label = "Strain"
		self.strain_index = 0
		self.isolates = {}
		self.group_handler(csvfile)
	def group_handler(self,csv_file):
		with open(csv_file, 'r') as csvfile:
			reader = csv.reader(csvfile)			
			for index, row in enumerate(reader):
				if index == 0:
					if self.group_label in row and self.strain_label in row:
						self.group_index = row.index(self.group_label)
						self.strain_index = row.index(self.strain_label)
					else:
						print("Group or Strain column not found in csv file. Exiting...")
						sys.exit(1)
				else:
					if row[self.group_index] in self.isolates:
						self.isolates[row[self.group_index]].append(clean_name(row[self.strain_index]))
					else:
						self.isolates[row[self.group_index]] = [clean_name(row[self.strain_index])]
				
	def dictonary_reader(self, yaml):
		isolates = []
		name = list(yaml)[0]
		items = list(yaml.values())[0]
		for group in items:
			if isinstance(group,dict):
				_isolates = self.dictonary_reader(group)
				isolates = isolates + _isolates
			else:
				group_isolate = self.isolates[group]				
				isolates = isolates + group_isolate
				self.dictonary[group] = group_isolate
		self.dictonary[name] = isolates
		return isolates
	def reader(self,lista):
		for item in lista:
			if isinstance(item,dict):
				self.dictonary_reader(item)
			else:
				self.dictonary[item] = self.isolates[item]
		return self.dictonary



