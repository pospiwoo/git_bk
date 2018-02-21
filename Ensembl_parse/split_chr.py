import sys, os
with open(sys.argv[1],'r') as inFile:
	for line in inFile:
		if line.startswith('>'):
			oFile = open(line[1:].strip().replace(' ','')+'.fa','w')
			oFile.write(line)
			continue
#		print line
		oFile.write(line)

