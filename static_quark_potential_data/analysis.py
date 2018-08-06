import os 


os.system("ls L16_Lt20_beta2.2977*.bin > dummy.txt") 
fin = open("dummy.txt")

inlist = fin.readlines()

fin.close()

for i in range(len(inlist)):
    inlist[i] = (inlist[i])[:-1]

outstr = "./autocorr 3"
for j in range(len(inlist)):
    outstr += " "
    outstr += inlist[j]
    
os.system(outstr)
os.remove("dummy.txt")

