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

values = list()

for j in range(len(inlist)):
    outstr = "Rscript fit_plot.R "
#    print(j)
#    outstr += str(11 + j)
#    outstr += ".bin "
    outstr += inlist[j]
    outstr += " > fit_param.txt" 
#    outstr += inlist[2*j+1]
#    print(outstr)
    print(inlist[j])
    os.system(outstr)
    smove = "mv Rplots.pdf "
    smove += inlist[j]
    smove += ".pdf"
    os.system(smove)
    fitin = open("fit_param.txt")
    outlist = fitin.readlines()
    values.append([float(inlist[j][1:-4]), float(outlist[1].split()[2]), float(outlist[2].split()[2])])
    sremove = "rm fit_param.txt"
    os.system(sremove)
values.sort()
data_file = open("potential.csv", "w")
for j in range(len(values)):
    data_file.write(str(values[j][0]))
    data_file.write(" ")
    data_file.write(str(values[j][1]))
    data_file.write(" ")
    data_file.write(str(values[j][2]))
    data_file.write("\n")


#print(values)  
os.remove("dummy.txt")

