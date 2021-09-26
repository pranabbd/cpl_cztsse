############# A python code to generate chemical potential landscapes (CPLs) for CuZnSn(S_1-xSe_x)4 ##################
############# Written by PRANAB SARKER ############# Email: srpranab@gmail.com #############

from   projCPL import CPL       #available at 
import matplotlib.pyplot as plt

#set variables for three alloy compositions of CZTSSe (primary phase)
x0_375   = 'Cu2Zn1Sn1S2.5Se1.5'
x0_5     = 'Cu2Zn1Sn1S2Se2'
x0_625   = 'Cu2Zn1Sn1S1.5Se2.5'
########################## USER INPUTS #####################################
inputData = {
    'pPhase' : [x0_375, x0_5, x0_625],      #CZTSSe with three different alloy composition 
    #formation energies (dHf) for CZTSSe phases
    #'pdHf'   : [-4.301, -4.227, -4.143],   #pristine CZTSSe
    'pdHf'   : [-4.169, -4.096, -4.022],    #Cu-vacant CZTSSe
    #'pdHf'   : [-4.276, -4.201, -4.117],   #(Cu_Zn + Zn_Cu)-CZTSSe
    #'pdHf'   : [-4.334, -4.263, -4.183],   #Zn_Cu-CZTSSe
    #'pdHf'   : [-4.237, -4.161, -4.084],   #Zn_Sn-CZTSSe
    #formation energies (dHf) for secondary phases
    'sPhase' : ['Cu2Zn1Sn1S4','Cu2Sn1S3', 'Cu2S1', 'Cu1S1', 'Zn1S1', 'Sn1S1', 'Sn1S2', 'Sn2S3', 'Cu2Zn1Sn1Se4', 'Cu2Sn1Se3', 'Cu2Se1', 'Cu1Se1', 'Zn1Se1', 'Sn1Se1', 'Sn1Se2'],
    'sdHf'   : [-4.587, -2.58, -0.807, -0.513, -1.895, -1.048, -1.328, -2.416, -3.935, -2.245, -0.558, -0.395, -1.621, -1.007, -0.827],
    #chemical species
    'sp'     : ['Cu', 'S', 'Se', 'Sn', 'Zn'],
    #chemical potential reference systems (Cu, S, Se, Sn, Zn)
    'mu'     : [-1.322, -4.1, -3.491, -3.835, -0.573], 
    #choose axes for CPL analysis
    'axes'   : ['Zn', 'Sn', 'Se'] # axis-order: x, y, diag 
}
######################### SET CONSTRAINTS ################################ 
fixed_dmu = {
    'sp'     : ['Cu', 'S', 'Se'], #change chemical species according to your need
    'dmu'    : [-1.06, -0.172, 0] #change these values 0 <= del_mu_sp <= pdHf (FE for CZTSSe) to set new constranints
    #'dmu'    : [-1.06, 0.0, -1.11]
    #'dmu'    : [-0.98, 0.0, -0.431]
}
######################### CALCULATE SECONDARY PHASE-LINE PASSING POINTS ON CPL ################################ 
#set the CPL object 
cpl = CPL(inputData, fixed_dmu)
#set constraints or not: 'yes' = constraints; 'none' = no constraints
#get all the primary phase data
pData = cpl.set_constraints('pPhase', 'pdHf', 'yes')
#get all corresponding secondary phases data
sData = cpl.set_constraints('sPhase', 'sdHf', 'yes')
#set axes: x-axis followed by y-axis
#set axes variables
xaxis = inputData['axes'][0] #Zn
yaxis = inputData['axes'][1] #Sn
diag  = inputData['axes'][2] #diagonal (anion-axis)
axisData = cpl.set_axes_coeff_lists(pData,xaxis, yaxis, diag)
#set all matrices of the same size to solve linear thermodynamic (equilibrium) equations for all phases
leData = cpl.set_linEq_data_lists(pData, sData, axisData)
#find intersecting points on X and Y axes
resultsXY = cpl.solve_with_axes(leData, axisData, xaxis, yaxis)
#find intersecting points on diagonal 
resultsD = cpl.solve_with_diag(leData, axisData, diag)
#write all intersecting points to the 'cpl.dat' file
cpl.write_cpl_data(resultsXY, resultsD)
######################### PLOT CPL ################################
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(1,3)
axes    = plt.gca()
for p  in range(len(resultsXY['pPhase'])):
    ax1   = ax[p].twiny().twinx()
    xorig = axisData['axceff'][p][0][1]
    yorig = axisData['axceff'][p][1][1]
    xmin  = axisData['axceff'][p][2][1]
    ymin  = axisData['axceff'][p][2][1]
    orig  = [xorig, yorig]
    xax   = [xmin, xorig]
    yax   = [yorig, ymin]
    diag1 = [xmin, yorig]
    diag2 = [xorig, ymin]
    sCounted = []
    #set axes limits
    #primary-phase x,y
    ax[p].set_xlim((xmin, xorig))
    ax[p].set_ylim((ymin, yorig))
    #secondary-phase x,y
    ax1.set_xlim((xmin, xorig))
    ax1.set_ylim((plot axisymin, yorig))
    #set axes lables
    ax[p].set_xlabel(r'$\mu_{Zn}$ (eV)', fontsize=22)
    ax[p].set_ylabel(r'$\mu_{Sn}$ (eV)', fontsize=22)
    ax[p].set_title(label=resultsXY['pPhase'][p], fontsize=22, color="blue")
    #customize tick
    ax[p].tick_params(direction='out', axis='both', length=6, width=1.5, colors='black',
               grid_color='r', grid_alpha=0.5)
    ax1.tick_params(direction='out', axis='both', length=6, width=1.5, colors='black',
               grid_color='r', grid_alpha=0.5)
    axes.tick_params(direction='out', axis='x', length=6, width=1.5, colors='black',
               grid_color='r', grid_alpha=0.5)
    #set axis_size
    for side in ax[p].spines.keys():  # 'top', 'bottom', 'left', 'right'
        ax[p].spines[side].set_linewidth(1.5)
    #plot axes
    ax1.plot(orig, xax, color='black')
    ax1.plot(yax, orig, color='black')
    ax1.plot(diag1, diag2, color='black')
    for s  in range(len(resultsXY['sPhase'][p])):
        #plot secondary-phase lines
        if s < (len(resultsXY['sPhase'][p]) - 1) and resultsXY['sPhase'][p][s] == resultsXY['sPhase'][p][s+1]:
            sPhaseXY = resultsXY['sPhase'][p][s]
            sCounted.append(sPhaseXY) 
            xval1  = round(resultsXY['points'][p][s][0], 3)
            xval2  = round(resultsXY['points'][p][s+1][0], 3)
            yval1  = round(resultsXY['points'][p][s][1], 3)
            yval2  = round(resultsXY['points'][p][s+1][1], 3)
            pt1    = [xval1, xval2]
            pt2    = [yval1, yval2]
            ax1.plot(pt1, pt2, label = sPhaseXY)
            ax1.legend(loc=3, prop={'size': 10})
        #plot secondary-phase lines
        else:
            sPhase = resultsXY['sPhase'][p][s]
            if sPhase not in sCounted: 
                xval1  = round(resultsXY['points'][p][s][0], 3)
                xval2  = round(resultsD['points'][p][s][0], 3)
                yval1  = round(resultsXY['points'][p][s][1], 3)
                yval2  = round(resultsD['points'][p][s][1], 3)
                #print(sPhase, xval1, xval2)
                #print(sPhase, yval1, yval2)
                pt1    = [xval1, xval2]
                pt2    = [yval1, yval2]
                ax1.plot(pt1, pt2, label = sPhase)
                ax1.legend(loc=3, prop={'size': 10})

plt.show()
