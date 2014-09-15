# file: $OPENFOAM_CASE_DIR/postProcessing/plot_forceCoeffs.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: plots aerodynamic coefficients of the flying snake

# Plot aerodynamic coefficients of the flying snake
# Olivier Mesnard
# mesnardo@gwu.edu

import matplotlib.pyplot as plt
import sys
import os

def getMeanValue(t,val):
	### To change if needed ###
	tStart,tEnd = 20.0,80.0
	### - ###
	mean = 0.
	index = 0
	for i in range(len(t)):
		if (tStart<=t[i]<=tEnd):
			mean += val[i]
			index += 1
	if (index != 0): mean = mean/index
	else: mean = None
	return mean

def main(arg):

	### To change if needed ###
	Re = 1000
	AOA = 25
	meshSize = '450K'
	### - ###

	pwd = os.path.dirname(os.path.abspath(__file__))

	fileList = sorted(os.listdir(pwd+'/forceCoeffs'))
	t,Cm,Cd,Cl = [],[],[],[]
	for i in range(len(fileList)):
		inFile = open(pwd+'/forceCoeffs/'+fileList[i]+'/forceCoeffs.dat')
		for line in inFile:
			data = line.split()
			if ('#' in data):
				pass
			else:
				t.append(float(data[0]))
				Cm.append(float(data[1]))
				Cd.append(float(data[2]))
				Cl.append(float(data[3]))
		inFile.close()

	Cd_mean = getMeanValue(t,Cd)
	Cl_mean = getMeanValue(t,Cl)
	print 'OpenFOAM - mean values: Cd=',Cd_mean,' - Cl=',Cl_mean

	plt.figure(0)
	plt.grid(True)
	plt.xlabel(r'time (s)',fontsize=12)
	plt.ylabel(r'Coefficients',fontsize=12)
	myLegend = []
	plt.plot(t,Cd,'b-',linewidth=2)
	myLegend.append(r'$C_d$ - icoFoam '+meshSize)
	plt.plot(t,Cl,'g-',linewidth=2)
	myLegend.append(r'$C_l$ - icoFoam '+meshSize)

	if ('--cuIBM' in arg):
		t_cuIBM,Cd_cuIBM,Cl_cuIBM = [],[],[]
		inFile = open(pwd+'/forceCoeffs_Re'+str(Re)+'_AOA'+str(AOA)+'_cuIBM','r')
		for line in inFile:
			data = line.split()
			t_cuIBM.append(float(data[0]))
			Cd_cuIBM.append(2.*float(data[1]))
			Cl_cuIBM.append(2.*float(data[2]))
		inFile.close()
		plt.plot(t_cuIBM,Cd_cuIBM,'k-.',linewidth=2)
		myLegend.append(r'$C_d$ - cuIBM')
		plt.plot(t_cuIBM,Cl_cuIBM,'k--',linewidth=2)
		myLegend.append(r'$C_l$ - cuIBM')
		Cd_mean_cuIBM = getMeanValue(t_cuIBM,Cd_cuIBM)
		Cl_mean_cuIBM = getMeanValue(t_cuIBM,Cl_cuIBM)
		print 'cuIBM - mean values: Cd=',Cd_mean_cuIBM,' - Cl=',Cl_mean_cuIBM
	plt.title('--- Flying Snake - Re='+str(Re)+' - AOA='+str(AOA)+' ---')
	plt.legend(myLegend,'best',prop={'size':12})
	plt.annotate(r'<$C_l$>='+str(Cl_mean)+'\n<$C_d$>='+str(Cd_mean),xy=(0,0.6),xytext=(0,0.6),\
				bbox=dict(fc='white'))
	plt.xlim(t[0],t[-1])
	plt.ylim(0.5,3.)
	if ('--save' in arg):
		plt.savefig(pwd+'/forceCoeffs_'+str(Re)+'_'+str(AOA)+'.png')
	if ('--show' in arg): plt.show()
	
#---------------------------------------

#---------------------------------------
if (__name__ == '__main__'):
	print '\nplotting aerodynamic coefficients...\n'
	main(sys.argv)
#---------------------------------------
