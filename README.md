# cpl_cztsse
A python code to generate chemical potential landscapes for CZTSSe.
NOTE: make sure calcCPL.py, projCPL.py, and run_cpl.sh are in the same directory

How to run:
1. install python 3+
2. install numpy, scipy, and matplotlib (for plotting) packages via the following commands
   pip install numpy (or 'conda install numpy' if anaconda is installed on your computer)
   pip install scipy (or 'conda install scipy' if anaconda is installed on your computer)
   pip install matplotlib (or 'conda install matplotlib' if anaconda is installed on your computer)
3. type './run_cpl.sh' and hit ENTER OR type 'python calCPL.py > out.dat' and hit ENTER
4. if the matplotlib is used, please note that it is for quick visulization only; adjust it for better visulization
5. open 'out.dat' to see which secondary phases are not present in the CPL at a given thermodynamic conditions
6. open 'cpl.dat' file to see all chem. pot. points for all phases present in the CPL at a given thermodynamic conditions


How to modify constraints:
1. open calCPL.py file
2. go to the follownig section

######################### SET CONSTRAINTS ################################
fixed_dmu = {
    'sp'     : ['Cu', 'S', 'Se'], #change chemical species according to your need
    'dmu'    : [-1.06, -0.172, 0] #change these values within 0 <= del_mu_sp <= pdHf (CZTSSe formation energy) to set new constranints
}

3. change chem. species ('SP') or chem. pot. ('dmu') from the dictionary lists
4. make sure the order  of chem. species name(s) and the correponding value(s) is same
