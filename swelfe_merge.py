#!/usr/bin/python
import string, re, sys, os 


if len(sys.argv) <4 :
	sys.exit ("utilistation : ./supprime_chevauch.py fichier.swelfe pourcentage fichier.sortie ")

fich = open(sys.argv[1] , 'r')
lignes = fich.readlines()
fich.close()

pourc = float(sys.argv[2])

fich = open(sys.argv[3] , 'w')

ecrit = 1
ancien_nom = ""
for i in lignes  :  
	if i[0:2] == ' @':
		fich.write( i,)
	elif i[0] == '#' or i[0] == '=':
		if ecrit == 1 :
			fich.write( i,)
	elif i[0] == '\n':
		fich.write( i,)
	elif i[0] == '@' :
		if ecrit == 1 :
			fich.write( i,)
	elif i[0] == '*' :
		fich.write( i,)

	else : 
		tmp = string.split(i)
		nom = tmp[0] 
		if nom != ancien_nom :
			seq1 = set(range(int(tmp[3]), int(tmp[3]) + int(tmp[10])))
			seq2 = set(range(int(tmp[4]), int(tmp[4]) + int(tmp[10])))
			ancien_nom = nom
			tab1 = []
			tab1.append(seq1)
			tab2 = []
			tab2.append(seq2)
			tab = []
			tab.append(tmp)
			fich.write( i,)
			ecrit = 1
		else : 
				
			seq11 = set(range(int(tmp[3]), int(tmp[3]) + int(tmp[10])))
			seq22 = set(range(int(tmp[4]), int(tmp[4]) + int(tmp[10])))
			flag = 0
			for j in range(len(tab1)) :
				if int(tmp[10])  > int(tab[j][10]) :
					lg = int(tab[j][10])
				else :
					lg =int(tmp[10])
				if len(tab1[j] & seq11) >= pourc * lg  and len(tab2[j] & seq22) >= pourc * lg : 
					flag = 1
					ecrit = 0
			if flag == 0 :				
				tab1.append(seq11)
				tab2.append(seq22)
				tab.append(tmp)
				fich.write( i,)
				ecrit = 1
	


fich.close()		
		 
