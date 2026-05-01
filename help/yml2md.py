#!/usr/bin/env python3
import yaml
import io
import sys

with io.open(sys.argv[1], 'r') as stream:
        data = yaml.safe_load(stream)

sections = ["name","title","description","details","authors","see_also","example","references"]

for entry in data:
	with io.open("md/"+entry['name']+".md", 'w') as outfile:
		for section in sections:
			outfile.write("## "+section+"\n")
			if entry[section] is not None:
				if section == 'example':
					entry[section] = "\t" + entry[section].replace('\n', "\n\t")
				
				if type(entry[section]) == list:
					for item in entry[section]:
						if type(item) == dict:
							yaml.dump(entry[section], outfile)
							break
						else:
							outfile.write(item+"\n")
				else:
					outfile.write(entry[section]+"\n")
