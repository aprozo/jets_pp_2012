
# Jets in pp 2012 200 GeV

This is an analysis for work in special container which analyzes `TStarJetPicoDsts` https://github.com/wsu-yale-rhig/TStarJetPicoMaker.

The companion external repository with how to read it
https://github.com/kkauder/eventStructuredAu 


The image container is `star_star.simg`.
- You need to install either [Docker engine](https://docs.docker.com/get-started/get-docker/) or [Apptainer (singularity)](https://apptainer.org/docs/admin/main/installation.html).
For simplier Apptainer (singularity) installation:
```bash
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

### Notes
File `/gpfs01/star/pwg/youqi/run12/embedding/P12id/picos/20235003/out/pt-hat1115_000.root` does not contain Mc entries




## Instruction on how to produce TStarJetPicoDsts from minimcs and MuDsts:

* 1. Make a List of minimcs and MuDsts
**Command for making the full list:**
    ```bash
    find "$PWD" -type f | sort >> ~/TStarJetPicoMaker/<name_of_list>.list
    ```
**Note:** The `| sort` flag is necessary even if your inputs seem ordered by eye, because `find` will traverse the directory tree in the order items are stored within the directory entries. This will (mostly) be consistent from run to run on the same machine and will essentially be "file/directory creation order" if there have been no deletes.  
However, some file systems will re-order directory entries as part of compaction operations or when the size of the entry needs to be expanded, so there's always a small chance the "raw" order will change over time. If you want a consistent order, feed the output through an extra sorting stage.
    
**Command for separating MuDsts and minimcs:**
    ```bash
    sed '/minimc/!d' <name_of_list>.list >> <minimc_list_name>.list
    ```
    And similar for MuDsts.

* 2. Split Each List into Equal Files/Lines
Use the `split` command. Example:
```bash
    split --suffix-length=2 --numeric-suffixes --lines=100 --additional-suffix=.list --verbose <input_filename> <output_prefix>_
```
This splits into files with 100 lines each.  
Or, replace `--lines=100` with, e.g., `-n l/5` to split evenly into 5 files.
    

* 3. Change Filename Base in `macros/MakeTStarJetPico_example.cxx`

(Around line 100 currently.)  
    If your lists have a different structure than e.g. `MuDsts1115_00` (where `1115` is the pT-hat range and `00` is the first list of 100 MuDsts), you may need to rewrite this block.

Test in an interpreter (e.g., a ROOT environment) until you get a sensible value for `std::string unique_name`.

* 4. Change the Call in `submit/jetPicoProduction_example.xml`

Match that number of lines. (Nominally, I have 100 files per job, so there will probably be a `100` as the penultimate argument to be changed to suit your job size.)
    

* 5. Adjust `find_mu`, `find_mc` and `slice_name` in `submit/submit.py`

Adjust to match the lists of files you're inputting, and update `slice_name` so the XML script knows the name of the files it's copying from `$SCRATCH` to local.
    
The way `slice_name` is coded, be careful if you have lists like `minimcs1115.list` before splitting into 100 lines each. The Python looks for `minimcs`, so name these instead `minimc1115.list` to avoid double-counting.


* 6. Adjust Default Arguments in `submit/submit.py`

Update to suit your case so the defaults will work and you won’t have to input them each time.
    

* 7. Submit Jobs

Once the above is adjusted, simply run:
    
```bash
    python submit/submit.py
```
    
This will use the XML script `n` times, where `n` is the number of input minimc or MuDsts lists.
    


* 8. If You Adjust Anything in `StRoot/TStarJetPicoMaker/`

Remember to run:
```bash
./macros/compile.csh
```


* 9. Ensure `log/tmplogs` Directory Exists

Make sure the directory exists where the submit script expects it, or it will try to submit jobs until the very end, then fail.
