# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 14:29:57 2016

@author: swoo
"""
import os
from os.path import dirname, abspath

class Path:
    def __init__(self):
            
        self.script_dir = dirname(__file__)
        self.parent_dir = dirname(dirname(dirname(abspath(__file__)))) # py2exe
        self.parent_dir = dirname(dirname(abspath(__file__))) # normal
        
        self.out_dir = os.path.join(self.parent_dir, 'Output')
#        if not os.path.exists(self.out_dir):
#            os.makedirs(self.out_dir)
            
        self.input_dir = os.path.join(self.parent_dir, 'Input')
        self.encoding_dir = os.path.join(self.parent_dir, 'encoding')
        
    def makeNewOutputDir(self, new_output_dir):
        self.out_dir = os.path.join(self.out_dir, new_output_dir)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

            
###################### Deprecated ######################
#        try:
#            self.script_dir = dirname(__file__)
#            self.parent_dir = dirname(dirname(dirname(abspath(__file__))))
#            afile = open(os.path.join(os.path.join(self.parent_dir, 'test.txt')), mode='w')
#            afile.write(self.script_dir)
#            afile.write("\n")
#            afile.write(self.parent_dir)
#            afile.close()
#        except NameError:  # We are the main py2exe script, not a module
#            self.script_dir = dirname(sys.argv[0])
#            self.parent_dir = dirname(dirname(abspath(sys.argv[0])))
#            afile = open(os.path.join(os.path.join(self.parent_dir, 'test.txt')), mode='w')
#            afile.write(self.script_dir)
#            afile.write("\n")
#            afile.write(self.parent_dir)
#            afile.close()
#            #self.parent_dir = os.path.pardir(abspath(sys.argv[0]))