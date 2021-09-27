# A python code to generate chemical potential landscapes for CZTSSe.

NOTE: make sure calcCPL.py, projCPL.py, and run_cpl.sh are in the same directory

## How to run:
- install python 3+
- install numpy, scipy, and matplotlib (for plotting) packages via the following commands, if these packages are not already installed
   `pip install numpy` (or `conda install numpy` if anaconda is installed on your computer)
   `pip install scipy` (or `conda install scipy` if anaconda is installed on your computer)
   `pip install matplotlib` (or `conda install matplotlib` if anaconda is installed on your computer)
- with the bash script: type `chmod a+x run_cpl.sh` and hit ENTER and then again type `./run_cpl.sh` and hit ENTER
- without the bash script: type `python calCPL.py > out.dat` and hit ENTER
- if the matplotlib is used, please note that it is for quick visulization only; customize it by clicking the 'Configure subplots' icon for optimal visulization
- open 'out.dat' to see which secondary phases are not present in the CPL at a given thermodynamic conditions
- open 'cpl.dat' file to see all chem. pot. points for all phases present in the CPL at a given thermodynamic conditions


## How to modify constraints:
- open calCPL.py file
- go to the follownig section
 
```python
######################### SET CONSTRAINTS ################################

fixed_dmu = {
    'sp'     : ['Cu', 'S', 'Se'], #change chemical species according to your need
   
   'dmu'    : [-1.06, -0.172, 0] #change these values within 0 <= dmu <= pdHf 
                                 #(pdHf = CZTSSe formation energy) to set new constranints

}
```

- change chem. species ('SP') or chem. pot. ('dmu') from the dictionary lists
- make sure the order  of chem. species name(s) and the correponding value(s) is same
