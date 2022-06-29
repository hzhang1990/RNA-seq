#This program genrates the length of a Cigar string in the last coloum of the ".sam" file.
import sys

file_o1 = open(sys.argv[2], 'w')  ## output file, expmple : accepted_hits_end.sam

file1 = open(sys.argv[1])  ## input file, exmpale : accepted_hits.sam
length = 0
dic = {}
for line in file1:
    if line != "\n":
        start = 0
        end = 0
        splits = line.split('\t')
        if splits[0].strip() not in ["@HD", "@PG", "@SQ"]:  ## to remove the header
            length = len(splits)
            start = int(splits[3].strip())
            sums = 0
            s = ''
            # to convert the cigar string into length (integer)
            for i in splits[5].strip():  # check for cigar string
                if i.isdigit():
                    s = s + i
                else:
                    sums = sums + int(s)
                    s = ''
            end = start + sums
            line1 = [line.strip(), '\t' + str(end), '\n']
            file_o1.writelines(line1)
file_o1.close()
