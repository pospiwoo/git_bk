# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 09:54:59 2016

@author: Sunghee Woo
"""
import os #, time, sys
import Tkinter, tkFileDialog #, Tkconstants, tkMessageBox
import Ultra_Short_Read_Pipeline as USR
import Path_modules as Path

class USR_pipeline_GUI(Tkinter.Frame):
    def __init__(self, root):
        self.path = Path.Path()
        self.fastq = '' # read FASTQ file: 'sample_2_high.fastq'
        self.target_fa = '' # target sequence FASTA file: 'target_sequeces.fa'
        self.output_dir = '' # output dir: ./Output/XXX        
        Tkinter.Frame.__init__(self, root)
        #self.file_opt = options = {} # Currently using default options

        # define buttons
        button_row_ind = 0
        Tkinter.Label(self, text='Select read FASTQ file').grid(row=button_row_ind,column=0)
        Tkinter.Button(self, text='Select read FASTQ file', command=self.ask_fastq_filename).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Label(self, text='Select target FASTA file').grid(row=button_row_ind,column=0)
        Tkinter.Button(self, text='Select target FASTA file', command=self.ask_target_seq_file_name).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Label(self, text='New output directory name').grid(row=button_row_ind,column=0)
        self.output_dir = Tkinter.StringVar()
        Tkinter.Entry(self, textvariable=self.output_dir).grid(row=button_row_ind,column=1)
        
        button_row_ind += 1
        Tkinter.Button(self, text ='Run', command = self.submit).grid(row=button_row_ind,column=0)

    def ask_fastq_filename(self):
        self.fastq = tkFileDialog.askopenfilename(title="Select read FASTQ file", initialdir=self.path.input_dir)
        Tkinter.Label(self, text=os.path.basename(self.fastq)).grid(row=0,column=2)

    def ask_target_seq_file_name(self):
        self.target_fa = tkFileDialog.askopenfilename(title="Select target FASTA file", initialdir=self.path.input_dir)
        Tkinter.Label(self, text=os.path.basename(self.target_fa)).grid(row=1,column=2)

    def submit(self):
        Tkinter.Label(self, text='Running').grid(row=3,column=0)
        self.path.makeNewOutputDir(str(self.output_dir.get()))
        USR_inst = USR.USRController(self.fastq, self.target_fa, self.path)
        USR_inst.Process()
        Tkinter.Label(self, text='Finished').grid(row=3,column=0)


if __name__=='__main__':
    root = Tkinter.Tk()
    USR_pipeline_GUI(root).pack()
    root.wm_title("Hyb&Seq pipeline GUI")
    root.mainloop()
    