import sys, os

max_start = 248945028+100
div_many = 10
start_term = int(max_start/div_many)

oFile_list = []
for i in range(0,div_many):
    oFile = open(sys.argv[1].replace('.gff3','_'+str(i)+'.gff3'),'w')
    oFile_list.append(oFile)


with open(sys.argv[1],'r') as inFile:
    for line in inFile:
        if line.startswith('#'):
            continue
        data = line.strip().split('\t')
        start = int(data[3])
        chr_str = data[0]
        if chr_str != '1':
            continue
        div =  int(start / start_term)
        oFile_list[div].write(line)

#        if start > max_start:
#            max_start = start
#print (max_start)


for i in range(0,div_many):
    oFile_list[i].close()



#1	ensembl_havana	CDS	234472714	234472811	.	-	2	ID=CDS:ENSP00000040877;Parent=transcript:ENST00000040877;protein_id=ENSP00000040877
# 248945028
