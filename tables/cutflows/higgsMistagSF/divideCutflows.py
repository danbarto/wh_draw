import sys
sys.path.insert(1, '/home/users/fechen/CMSSW_10_2_9/src/WH_studies/Tools/python/')
#sys.path.insert(1, '/home/users/fechen/CMSSW_10_2_9/src/WH_studies/')
#print sys.path
import u_float as uf
#from '/home/users/fechen/CMSSW_10_2_9/src/myWH_studies/Tools/python/' import u_float as uf
#from WH_studies.Tools.u_float import u_float as uf
import fileinput

# Changes False to True and True to False
def change(var):
    if var:
        return False
    else:
        return True

# Converts latex to uf. Ex: "yield$\\pm$uncertainty" -> uf(yield, uncertainty)
def texToUf(tex):
#    print tex
    y = float(tex.split('$\\pm$')[0][1:])
    unc = float(tex.split('$\\pm$')[1][:-1])
    return uf.u_float(y, unc)

# Converts uf to latex. Ex: "12312+-321" -> 12312$\\pm$321
def ufToTex(myuf):
    return str(round(myuf.val, 2)) + "$\\pm$" + str(round(myuf.sigma, 2))

# Un-tokenize ie re-combine a bunch of list elements
def untokenize(arr, delim = ' '):
    mystr = ""
    for elemCount, elem in enumerate(arr):
        if len(arr) == elemCount+1:
            mystr += (elem + "\\\\")
        else:
            mystr += (elem + delim)
    return mystr

    
# Check arguments
if 4 != len(sys.argv):
    print "Three arguments only!"
    sys.exit()

# Open input and output files
numerFile = open(sys.argv[1], 'r')
denomFile = open(sys.argv[2], 'r')
#outFile = open("outFile.tex", 'w')

# Holds yields
numerYields = []
denomYields = []

# Extracts yields
numerLines = numerFile.readlines()
found = False
for line in numerLines:
    if "\hline" == line.strip():
        found = change(found)
    if found:
#        print line.strip().split('&')
        numerYields.append(line.strip().split('&'))

# Extracts yields
denomLines = denomFile.readlines()
found = False
for line in denomLines:
    if "\hline" == line.strip():
        found = change(found)
    if found:
#        print line.strip().split('&')
        denomYields.append(line.strip().split('&'))

# Get rid of first array, which is just "hline"
numerYields.pop(0)
denomYields.pop(0)
# Get rid of "\\\\" in the last element of every array
for i in range(len(numerYields)):
    numerYields[i][-1] = numerYields[i][-1][:-2]
    denomYields[i][-1] = denomYields[i][-1][:-2]

# Ratios
ratios = []

# Remakes tex lines after dividing
for i in range(len(numerYields)):
    ratio = []
    for elemCount in range(len(numerYields[i])):
        if 0 == elemCount:
            ratio.append(numerYields[i][elemCount]) # The cut description
        else:
            # division done here
            yAndUnc = ufToTex(texToUf(numerYields[i][elemCount])/texToUf(denomYields[i][elemCount]))            
            ratio.append(yAndUnc)
    ratios.append(ratio)
for i in range(len(ratios)):
    ratios[i] = untokenize(ratios[i], ' & ')

for ratio in ratios:
    print ratio
#
# Generates a divided tex file
from shutil import copyfile
copyfile(sys.argv[1], sys.argv[3])
#print "printing new"

#with fileinput.FileInput(sys.argv[3], 
#                         inplace = True) as f: 
#for line in f:
#    if "pm" in line:
#        print 
#    for line in f: 
#        if "search this line and change it\n" == line: 
#            print("changing the matched line with this line", 
#                  end ='\n') 
#        else: 
#            print(line, end ='') 
count = 0
#for line in fileinput.input([sys.argv[3]], inplace=True):
#    if "pm" in line:
#        print line.rstrip().replace(line, ratios[count]),
#        #        sys.stdout.write(line)
#        #        print line
##        print(ratios[count], end = '\n')
#        count+=1
##    else:
##        print(line, end = '')
##    if line.strip().startswith('initial_mass = '):
##        line = 'initial_mass = 123\n'
##    sys.stdout.write(line)

with open(sys.argv[1], 'r') as input_file:
    with open(sys.argv[3], 'w') as output_file:
        for line in input_file:
            if "pm" in line.strip():
                output_file.write("    " + ratios[count] + "\n")
                count+=1
        #            output_file.write('new line\n')
            else:
                output_file.write(line)
