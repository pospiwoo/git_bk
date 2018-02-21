from distutils.core import setup
import py2exe
import os


script_dir = os.path.dirname(__file__)
log_file_name = os.path.join(script_dir, 'Ultra_Short_Read_Pipeline_GUI.py')

#options = {'py2exe': {
#           'compressed':1,  
#           'bundle_files': 2, 
#           'dist_dir': "my/dist/dir"
#           'dll_excludes': ['w9xpopen.exe']
#           }}



setup(console=[log_file_name],
		options = {'py2exe': {
		'compressed':1,
	   'dist_dir': "C:\\Work\\Compiled\\py2exe_redirect",
	   'dll_excludes': ['w9xpopen.exe']
	   }}
	   )
