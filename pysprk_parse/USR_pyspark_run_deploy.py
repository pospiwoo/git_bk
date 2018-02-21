# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:21:35 2017

@author: Sunghee Woo
"""
import sys, os, time
import Assembly_pyspark_deploy_modules as Assm
import Path_modules as Path
#from random import random
#from operator import add
from pyspark import SparkContext
from datetime import timedelta



class USRPysparkController:
    def __init__(self, fastq, target_fa, path_inst):
        self.Assm_view = Assm.AssemblyView(path_inst)
        self.Assm_model = Assm.AssemblyModel(fastq, target_fa, self.Assm_view)
        self.Assm_controller = Assm.AssemblyController(self.Assm_model, self.Assm_view)
        
    def Process(self):
        start = time.time()
        self.Assm_controller.Process()
        end = time.time()
        self.Assm_view.print_time(start, end)
        
def print_time(start, end):
    time_in_sec = end - start
    print_str = 'Excecution time: ' + str(timedelta(seconds=time_in_sec))
    print print_str
    return
    
def parse_line(line_str):
    return_str = []
    return_str.append('1')
    return_str.append('2')
    return_str.append('3')
#    return_str.append('# ')
#    return_str += line_str.split(',')
    return return_str
#    return ['1', '2', '3']
#    start = time.time()
#    self.Assm_controller.Process()
#    end = time.time()
#    self.Assm_view.print_time(start, end)

if __name__=='__main__':
    path = Path.Path()
    if len(sys.argv) < 3:
        print('Insufficient arguments')
        print('Usage: python XXX.py <FASTQ> <Target_FASTA>')
        exit()
    elif len(sys.argv) >= 3:
        fastq = sys.argv[1]
        target_fa = sys.argv[2]
        
        Assm_controller = Assm.AssemblyPysparkProcessing(fastq, target_fa)
        sc_inst = SparkContext(appName="USR_spark")
        sc_inst.addPyFile('modules.zip')
        Assm_controller.Process(sc_inst)
        sc_inst.stop()
#        Assm_controller.mergeFASTA()
#        Assm_controller.mergeCOV()
#        Assm_controller.mergeFASTAandSumCov()
#        Assm_controller.writeInfo()
        
#        sc = SparkContext(appName="USR_spark")
#        text_file = sc.textFile(fastq)
#        file_data = text_file.map(parse_line)
##        file_data = text_file.flatMap(lambda line: line.split("\n"))
##        text_file = sc.textFile(file_data)
#        #file_data.take(15)
#        file_data.saveAsTextFile(target_fa)