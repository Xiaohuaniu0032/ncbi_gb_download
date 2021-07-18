from Bio import SeqIO
record = SeqIO.read("AB981583.gb", "genbank")

print(record.annotations)
for feature in record.features:
	#print(feature)
	if feature.type == "source":
		#strain = feature.qualifiers.get('strain','NA')
		#print(strain)
		print(feature)
