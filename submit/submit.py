from __future__ import print_function

from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
import os
import argparse
import re
import subprocess
import time
import math
import numpy as np

def main(args):

  ## check xml file exists
  xml_file = os.path.join(os.getcwd(), args.submitscript)
  if not os.path.isfile(xml_file) :
    print('xmlfile doesnt exist!')
    return

  ## create output from input variables
  out_base = args.outputroot
  log_dir = os.path.join(out_base, "log")
  out_dir = os.path.join(out_base, "out")

  ## create log and output directories
  if not os.path.exists(log_dir):
    os.makedirs(log_dir)
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)

  print('Submitting jobs - parameters')
  print('log directory: ', log_dir)
  print('output directory: ', out_dir)
  
  file_list = []
  
  with open('data_file.list', 'r') as data_file_list:
    for data_file in data_file_list:
	file_list.append(data_file.strip())

  for data_file in file_list:

    print("submitting job for data file: ", str(data_file))

    submit_args = 'out=' + out_dir
    submit_args = submit_args + ',log=' + log_dir
    submit_args = submit_args + ',data_file=' + str(data_file)
    data_name   = os.path.basename(data_file)
    submit_args = submit_args + ',data_name=' +  str(data_name)
   
    star_submit = 'star-submit-template '
    star_submit = star_submit + '-template ' + xml_file
    star_submit = star_submit + ' -entities ' + submit_args
    
    ret = subprocess.Popen(star_submit, shell=True)
    ret.wait()
    if ret.returncode != 0:
      print('warning: job submission failure')

    
    
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Submit jobs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--submitscript', default='submit/submit.xml', help='the xml file for star-submit-template')
  parser.add_argument('--outputroot', default='/gpfs01/star/pwg/prozorov/jets_pp_2012/output', help='root directory for all output and logs')
  args = parser.parse_args()
  main( args )

