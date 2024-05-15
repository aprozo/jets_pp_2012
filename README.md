<<<<<<< HEAD
anytime i type something in “$ “ it means that you type whatever is after the $ in your terminal. 

1. login to RCF 
2. get a working node - "rterm -i”
3. convert your shell to bash - “bash”
4. remove any older copies of my working directory /gpfs01/star/pwg/elayavalli/ppRun12_analysis_code 
5. Copy the directory into your working directory 
6. run the singularity container “$ source runimage.sh”
7. if this worked, you should see your name become Singularity> at the very start of the terminal  
8. setup the appropriate libraries and paths inside the container  “$ source setup.sh"
9. the main analysis codes are available inside the ‘src’ directory - i created a new analysis called pptest which you can see the corresponding 
10. first cleanup the existing binaries by running “$ make clean” from your ppRun12_analysis_code directory 
11. Then compile the executable - “$ make”
12. you should see that you have one executable inside ’bin/‘ now 
13. there are three .txt files which include the necessary commands you need to run the analysis and the corresponding switches on pythia, geant, data. this is, in my view, truly the advantage of such a code in that your analysis method is 100% common and everything is controlled by the switches 
14. read the code and understand what calls what, how can you add new observables into the ResultStructure etc… add your own observables and you can save it to a tree or save it directly to histograms. You can also see how you can get the charged particles within the jet inside the analysis class’s method that runs every event. 
15. Have fun!  :) and as always, feel free to ask me any questions. 





The directory I use is /gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/submit. Inside the directory there is submit.py and submit.xml

submit.py - We have a for loop that loops over all the files in the file list (/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/data_file.list). One job will be submitted for each iteration of the loop. In each loop, we defined a string "submit_args" which contains arguments out, log, data_file and they will be passed to submit.xml

submit.xml - Line 35 allows us to enter the container. This is the equivalent of "source runimage.sh" when we run the code interactively. The options -B are directories that we want to be visible inside the container. It calls container.sh with singularity and reads the arguments out, log, data_file. data_file will then be passed to container.sh. Line 36 moves the output files of the jobs from a temporary location to the out directory you specified in the python script. Starting from line 48 is sandbox. I included all the files needed to run the jobs here, even the 65 GB of data (copied from Raghav's directory to mine), which is not the best way of doing it but makes my life easier.

container.sh (/gpfs01/star/pwg/youqi/run12/ppRun12_analysis_code/container.sh) - this contains the command lines to run the code and takes in arguments from submit.xml. $1 corresponds to the data_file variable in submit.xml

To run the code, create a /submit directory in your work directory, copy submit.py and submit.xml files from my directory, edit the code in those files (so that they read the file list you have, save output to your directory, etc). Create directories /scheduler/report, /scheduler/csh and /scheduler/list under /submit (see lines 70-72 of the xml file). Copy over and modify container.sh (for me this file was put outside the submit directory). Finally, do "python submit/submit.py" to submit the jobs.
=======
# jets_pp_2012
>>>>>>> c2dbe548387d7857e3a68874cb869431aa93fc64
