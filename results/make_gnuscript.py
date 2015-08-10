#!/usr/bin/env python3.4
import glob
flist = glob.glob("*.dat")

f = open("scrpt","w")
f.write("set view map\nset xlabel \"x\"\nset ylabel \"y\"\nset cblabel \"Temp. [*C]\" offset -6.5,0 rotate by 90\nset palette defined (0 \"royalblue\",1 \"turquoise\",3 \"yellow\",5 \"red\")\nset cbrange [0:100]\n")
flist = [ i.split("_")[0] for i in flist]
flist.sort(key=lambda x: float(x.split("_")[0]))

#print(flist)
for i in range(len(flist)):
	pngname="-".join(flist[i].split(".",1))[:6]
	element="set terminal png\nset output \"picts/%05d"%i+""+".png\"\nset title \""+flist[i] +"s\"\nsplot \""+flist[i]+"_out.dat\" w pm3d t \"\"\n"
	f.write(element)


