import os 

os.system("ls r*.csv > dummy.txt") 
fin = open("dummy.txt")

inlist = fin.readlines()

fin.close()
#inlist = os.listdir(./L*)

#print(inlist)

for i in range(len(inlist)):
    inlist[i] = (inlist[i])[:-1]

#print(inlist)


for j in range(len(inlist)):
    outstr = "Rscript fit_plot.R "
#    print(j)
#    outstr += str(11 + j)
#    outstr += ".bin "
    outstr += inlist[j]
    outstr += " > fit_param.txt" 
#    outstr += inlist[2*j+1]
#    print(outstr)
    os.system(outstr)
    fitin = open("fit_param.txt")
    outlist = fitin.readlines()
    print(outlist)
    os.remove("fit_param.txt")

    
os.remove("dummy.txt")

