#!/usr/bin/python
'''This program is used to sort many lines use one colum.'''
from __future__ import division
import os
class SortedByColum:
    def __init__(self,file_in,file_out,colum = 1):
        self.file_in = file_in
        self.file_out = file_out
        self.colum = colum
    def sorted_by_colum(self):
        f_h = open(self.file_in)
        array = f_h.readlines()  # this makes 'array' a list ,every line id one element
        array = [x.split('\t') for x in array]  #this make 'array' a list contains every colum whic is also a list
        array = [[int(x[1]),x[0],int(x[2].replace('\n','')),int(x[2])/int(x[1])]
                for x in array] #this replaced the '\n',and modified the position
        #array = [[int(x[1]),x[0],x[2]] for x in array]
        f_h.close()
        
        f_o = open(self.file_out,'w')
        array.sort()
        array = [str(x[1]) + '\t' + str(x[0]) + '\t' + str(x[2]) + '\t' + 
                str(x[3]) + '\n' for x in array]
        f_o.writelines(array)
        f_o.close()
        
        
        
if __name__ == '__main__':
    cur_dir = os.getcwd()
    file_dir = cur_dir[:-6]
    file_in = file_dir + '1018' + os.sep + 'parse' + os.sep + 'length_protein_domain_h'
    file_out = file_dir + '1018' + os.sep + 'parse' + os.sep + 'length_protein_domain_h_sorted_again'
    
    sbc = SortedByColum(file_in,file_out)
    sbc.sorted_by_colum()
    print "Done!"
       
  
