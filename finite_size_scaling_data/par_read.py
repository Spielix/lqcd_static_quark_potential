import os 

os.system("ls L*.bin > dummy.txt") 
fin = open("dummy.txt")

inlist = fin.readlines()

fin.close()
#inlist = os.listdir(./L*)

#print(inlist)

for i in range(len(inlist)):
    inlist[i] = (inlist[i])[:-1]

#print(inlist)


for j in range(len(inlist)):
    outstr = inlist[j]
    print(outstr)
    outstr = "./par_read "
    outstr += inlist[j]
    os.system(outstr)

os.remove("dummy.txt")

