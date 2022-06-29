### In this step you have to remove the ends genrated in step 2. Otherwise it will change
#file format of .sam-file, and can not be used for read count. Infile needs to be previously
#defined within the script. Outfile needs to be defindes in output-file.
import sys

file1 = open(sys.argv[1])  ## Inputfile from previous step
for line in file1:
    if line != "\n":
        splits = line.split('\t')
        if splits[0].strip() in ["@HD", "@PG", "@SQ"]:
            print line.strip()
        else:
            l = splits[0].strip()
            for i in range(1, len(splits) - 1):
                l = l + '\t' + splits[i].strip()
            print l
