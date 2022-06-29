# filter for cases like 1) same read id mapping to different position with same start and different end (on same chromosome) accept shorter in length
#			2) same read id mapping to different position  with different end and same end( same chromosome) accept shorter in length
#			3) same read id mapping at same position on same chromosome-> choose any of it

# this program uses the CIGAR string and NH field of .SAM file and removes the above mentioned error to avoid the bias in read counting.
# input of this program is from output of toprintend.py
import sys

file_o1 = open(sys.argv[2], 'w')  # output file name like accepted_hits_nonrepeated.sam

file1 = open(sys.argv[1])  # input file name like path/accepted_hits_end.sam
sdic = {}
dice = {}
for line in file1:
    if line != "\n":
        splits = line.strip().split('\t')
        start = int(splits[3].strip())
        chrom = splits[2].strip()
        end = int(splits[len(splits) - 1].strip())
        #		end=int(splits[5].strip())
        id = splits[0].strip()
        k = (id, chrom)
        kk = k + (start,)
        if kk not in sdic:
            sdic[kk] = [line.strip(), end - start + 1, end]
        elif kk in sdic:
            if end - start + 1 < sdic[kk][1]:
                sdic[kk] = [line.strip(), end - start + 1, end]

dice = {}
for i in sdic:
    if (i[0], i[1], sdic[i][2]) not in dice:
        dice[i[0], i[1], sdic[i][2]] = sdic[i]
    elif (i[0], i[1], sdic[i][2]) in dice:
        if sdic[i][2] - i[2] + 1 < dice[i[0], i[1], sdic[i][2]][1]:
            dice[i[0], i[1], sdic[i][2]] = sdic[i]

key = dice.keys()
key.sort()
for j in key:
    ln = [dice[j][0], '\n']
    file_o1.writelines(ln)

file_o1.close()
