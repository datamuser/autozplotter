#AutoZplotter is an autozygosity mapping algorithm which takes individual VCF files as input and enables visualisation and manual identification of regions that have longer stretches of homozygosity than would be expected by chance
#Source Paper: Erzurumluoglu et al, 2015. Identifying Highly Penetrant Disease Causal Mutations Using Next Generation Sequencing: Guide to Whole Process. Biomed Res Int. Article ID: 923491
#AutoZplotter is designed to handle VCF v4.0 (or similar formats)
#Make sure Python v2.7.13 (or 2.7.12) is installed
#Make sure to load XMing or a similar X11 display server (and 'X11 Tunnelling' enabled if using PuTTY on a remote server)
#Any questions, please contact Dr. A Mesut Erzurumluoglu (epmmee@bristol.ac.uk) or Dr. Tom R Gaunt (tom.gaunt@bristol.ac.uk)
#Usage: python AutoZplotter.py - A pop-up should appear where you can browse for and select the VCF file you'd like to visualise

#!/usr/local/bin/python
import sys
import matplotlib.pyplot as plt

from Tkinter import *
from tkMessageBox import *
from tkColorChooser import askcolor              
from tkFileDialog   import askopenfilename      

def calchet(winhet):
	sum = 0.0
	count = 0.0
	for item in winhet:
		sum += float(item)
		count += 1.0
	return sum/count
#try:
#regnfilename = sys.argv[1]
infilename = askopenfilename(defaultextension='.vcf',filetypes=[('Variant Call Format','*.vcf'),('All files','*.*')])#sys.argv[1]

#regnfile = open(regnfilename, 'r')
infile = open(infilename, 'r')
while 1:
	dataline = infile.readline()
	dataline = dataline.strip().split()
	if dataline[0] == "#CHROM":
		break
#results = []
#regions = regnfile.readlines()
#regnfile.close()

hzplotx = []
hzploty = []
hetplotx = []
hetploty = []

plotx = []
ploty = []

currentchr = "null"
counter = 0
windowhet = []

for i in range(25):
	plotx.append([])
	ploty.append([])

while 1:
	dataline = infile.readline()
	#print dataline
	dataline = dataline.strip().split()
	if len(dataline) < 5:
		break
	genocodes = dataline[9].split(":")
	chrom = dataline[0]
	if chrom in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","X","Y",1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]:
		chromn = chrom
	else:
		chromn = chrom[3:]
	if chromn == "X" or chrom == "X":
		chromn = 23
	elif chromn == "Y" or chrom == "Y":
		chromn = 24
	else:
		chromn = int(chromn)
	if chromn != currentchr:
		currentchr = chromn
		windowhet = []
		counter = 0
		
		
	pos = dataline[1]
	genotype = genocodes[0]
	if genotype == "0/0" or genotype =="1/1":
		thishet = 0.0
	else:
		thishet = 1.0
	if dataline[2][0:2] == "rs":
		windowhet.append(thishet)
	if counter > 20 and dataline[2][0:2] == "rs":
		null = windowhet.pop(0)
		plotx[chromn].append(pos)
		thisval = float(calchet(windowhet))
		ploty[chromn].append(float(chromn)*2+thisval)
	counter += 1
	
	readdepth = genocodes[1]
	#DP = genocodes[2]
	#GQ = genocodes[3]
	#PL = genocodes[4]
	#labelled = 0
	#for region in regions:
	#	region = region.strip().split()
	#	if region[1] == chrom and pos >= region[2] and pos <= region[3]:
	#		results.append([chrom,pos,dataline[2]+"("+dataline[3]+"/"+dataline[4]+")",genotype,readdepth,DP,GQ,PL,region[0],region[1],region[2],region[3]])
	#		labelled = 1
	#if labelled == 0:
	#	results.append([chrom,pos,dataline[2]+"("+dataline[3]+"/"+dataline[4]+")",genotype,readdepth,DP,GQ,PL,"NA","NA","NA","NA"])
	if (genotype == "0/0" or genotype =="1/1" or genotype == "0|0" or genotype =="1|1") and dataline[2][0:2] == "rs":
		hzplotx.append(pos)
		hzploty.append(float(chromn)*2+1+0.4)
	elif dataline[2][0:2] == "rs":
		hetplotx.append(pos)
		hetploty.append(float(chromn)*2+1+0.6)
	#print pos,chrom,genotype

		
		
		
		
fig = plt.figure()
ax = fig.add_subplot(111)

filelist = infilename.split("/")
print filelist
filetitle = filelist[-1]
fig.suptitle("File: " + filetitle, fontsize=11)

#for i in range(1,25):
#	ax.plot([0,250000000],[i,i],'k:')

plt.plot(hzplotx,hzploty,'ro',hetplotx,hetploty,'go')
#print "ready",plotx[0][0:5],ploty[0][0:5]

#ax.set_yticks((1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5))
ytick = []
for i in range(1,25):
	ytick.append(i*2+1)
ax.set_yticks(ytick)
ticklabels = []
for i in range(1,23):
	ticklabels.append("Chr "+str(i))
ticklabels.append("Chr X")
ticklabels.append("Chr Y")

labels = ax.set_yticklabels(ticklabels)

for i in range(len(plotx)):
	ax.plot([0,250000000],[(i+1)*2+1,(i+1)*2+1],'k:')
	ax.plot([0,250000000],[(i+1)*2,(i+1)*2],'k-')
	ax.plot(plotx[i],ploty[i],'b-')
print "ready",plotx[0][0:5],ploty[0][0:5]

#ax.axis('tight')
#plt.plot([1,4,8],[2,2,1])
plt.show()
print "done"
infile.close()
#for item in results:
#	print "\t".join(item)
		


#except:
#	print "\nscanautoz: A tool to scan autozygous regions from a vcf file\nUsage:\n\tscanautoz <vcf-file>\n"

