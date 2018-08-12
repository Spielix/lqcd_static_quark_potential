import os 

os.system("ls L12*.bin > dummy.txt") 
fin = open("dummy.txt")

inlist = fin.readlines()

fin.close()
#inlist = os.listdir(./L*)

#print(inlist)

for i in range(len(inlist)):
    inlist[i] = (inlist[i])[:-1]

#print(inlist)

outstr = "./autocorr_polya"

for j in range(len(inlist)):
    outstr += " 1 "
#    print(j)
#    outstr += str(11 + j)
#    outstr += ".bin "
    outstr += inlist[j]
#    outstr += " " 
#    outstr += inlist[2*j+1]
#    print(outstr)
    
os.system(outstr)
os.remove("dummy.txt")

