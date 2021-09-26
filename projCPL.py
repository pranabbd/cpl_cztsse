import collections
import itertools
import numpy as np
from   numpy import ones,vstack
from   scipy import linalg

class CPL(object):

    def __init__(self, inputs, constr):
        self.inputs  = inputs
        self.constr  = constr       #fixed chemical potential               
        self.cplData = collections.defaultdict(list)

    #plug the constrained value(s) in sumOf(delmu_i) = delHf equations 
    #get constrained sumOf(delmu_i) = delHf equations
    def set_constraints(self, whichPhase, whichFE, fixedVal):
        self.whichPhase = whichPhase 
        self.whichFE    = whichFE
        self.fixedVal   = fixedVal
        self.phaseData  = collections.defaultdict(list)
        for p in range(len(self.inputs[self.whichPhase])):
            stoic   = self.inputs[self.whichPhase][p] 
            #extract species names and stoichiometries for each phase
            items   = ["".join(x) for _, x in itertools.groupby(stoic, key=str.isalpha)]
            val     = 0
            ceff    = []
            var     = []
            #range(len(self.inputs[self.whichPhase]))
            for i in range(len(items)):
                #check whether or not the input has the correct format of stoichiometry: a symbol followed by a numeric
                if len(items[i]) > 2 and not items[i].isdigit:  #no numeric, only symbol
                    print('coefficient(s) missing for', stoic + str('.'), 'Check your inputs')
                    print('You must specify a numeric value e.g., AB = A1B1' + str('.'))
                    exit(0)
                #if there is no constraints
                if self.fixedVal == 'none':
                    val = 0.0
                    if items[i].isalpha():    #species strings
                        var.append(items[i])
                    else:
                        ceff.append(float(items[i])) #numeric strings (stoichiometries)
                #get total constraints
                elif self.fixedVal == 'yes':
                    if items[i] in self.constr['sp'] and items[i] in self.inputs['axes']: #e.g., Se is fixed and also an axis
                        idx = self.constr['sp'].index(items[i])
                        var.append(items[i])
                        ceff.append(float(items[i+1]))
                        val += float(items[i+1]) * self.constr['dmu'][idx]    #stoichiometry*delmu_i 
                    elif items[i] in self.constr['sp'] and items[i] not in self.inputs['axes']: #e.g., fixed-Cu but not an axis
                        idx = self.constr['sp'].index(items[i])
                        val += float(items[i+1]) * self.constr['dmu'][idx] 
                    #get variables and stoichiometries (coefficients)
                    elif items[i] not in self.constr['sp'] and items[i].isalpha():       # e.g,, Zn,Sn
                        var.append(items[i])
                        ceff.append(float(items[i+1]))
            pIdx   = self.inputs[self.whichPhase].index(stoic) 
            FE     = self.inputs[self.whichFE][pIdx] - val          #delHf - fixed-sumOf(delmu_i)
            #updated all constrained data into the dictionary
            self.phaseData['phase'].append(stoic)
            self.phaseData['formEn'].append(FE)
            self.phaseData['var'].append(var)
            self.phaseData['coeff'].append(ceff)

        return(self.phaseData)

    #set CPL axes/diagonal into the matrix form
    def set_axes_coeff_lists(self, pData, xaxis, yaxis, diag):
        self.xaxis  = xaxis 
        self.yaxis  = yaxis 
        self.diag   = diag 
        self.axData = collections.defaultdict(list)
        em          = []
        for p in range(len(self.inputs['pPhase'])):
            #set x-axis limit
            if self.xaxis in pData['var'][p]:
                idx  = pData['var'][p].index(self.xaxis)
                xmax = [pData['formEn'][p] / pData['coeff'][p][idx], 0]
            #set y-axis limit
            if self.yaxis in pData['var'][p]:
                idx  = pData['var'][p].index(self.yaxis)
                ymax = [0, pData['formEn'][p] / pData['coeff'][p][idx]]
            #set diagonal limit
            pdiag    = [xmax, ymax]
            #get slope and intercept for diagonal
            x, y     = zip(*pdiag)
            A        = vstack([x, ones(len(x))]).T       #y = mx + c; rewrite: y = Ap; A = [[x 1]]; p = [[m], [c]]
            m, c     = np.linalg.lstsq(A, y, rcond=None)[0]
            #get matrix for diagonal
            dxyCeff  = [-m, 1, 0]     # m     = coefficient of x; 1 for y; matrix for diagonal: [x, y] = [c]
            diaCeff  = [dxyCeff, c]   # list[0] = co-efficients of x and y; list[1] = c
            #get matrices for x-axis and y-axis
            xaxCeff  = [[0, 1, 0], 0]   # x-axis: 1y = 0
            yaxCeff  = [[1, 0, 0], 0]   # y-axis: 1x = 0
            self.axData['axis'].append([self.xaxis, self.yaxis, self.diag])
            self.axData['axceff'].append([xaxCeff, yaxCeff, diaCeff])
        
        return(self.axData)
    
    #transform all constrained FE equations into 3x3 matrices: [x-axis, y-axis, diagonal] = [FE]
    def set_linEq_data_lists(self, pData, sData, axisData):
        self.leData = collections.defaultdict(list)
        for p in range(len(self.inputs['pPhase'])):
            self.leData['pPhase'].append(self.inputs['pPhase'][p])
            self.leData['pSP'].append(pData['var'][p])
            self.leData['pceff'].append(pData['coeff'][p])
            self.leData['pFE'].append(pData['formEn'][p])
            pceff_tot = len(pData['coeff'][p])
            sec      =  []   #list for secondary phases
            ssp      =  []   #list for atom in secondary phases
            scoeff   =  []   #list for linear coefficients in secondary phases
            sfe      =  []   #list for formation energies for secondary phases
            for s in range(len(self.inputs['sPhase'])):
                sceff_tot = len(sData['coeff'][s])
                #declare a list according size of primary-phase
                #avoid phases whose all chemical species are constrained e.g., CuS when both chem. of Cu and S are fixed.
                scfs  = np.zeros(pceff_tot).tolist()  
                satms = np.zeros(pceff_tot).tolist()
                #avoid sce phases with  FE_sec < FE_primary; i.e., those sec phases will not sustain at the given thermodynamic conditions
                #no unconstrained chem. spices
                if sceff_tot == 0:# or sData['formEn'][s] >= 0:      
                    stoic = self.inputs['sPhase'][s]     
                    items   = ["".join(x) for _, x in itertools.groupby(stoic, key=str.isalpha)]
                    val = 0 
                    for i in range(len(items)):
                        if not items[i].isdigit:  #no numeric, only symbol
                            idx = self.constr['sp'].index(sData[items[i]])
                            val += float(items[i+1]) * self.constr['dmu'][idx] 
                    if val < pData['formEn'][p] or val >= 0:      
                        print('NOTE:')  
                        print('########################################################################################################################')  
                        print('Very poor (sumO_del_mu < del_Hf_primary) or non-equilibrium (sumOf_del_mu > 0) conditions for the following phase.')  
                        print( self.inputs['pPhase'][p],": ",'Sec. phase: ', self.inputs['sPhase'][s],"; ", 'FE_sec (eV): ', round(sData['formEn'][s], 3),\
                             '; FE_prim (eV): ', round(pData['formEn'][p], 3)) #sData['var'][s])     
                        print('The above phase will not co-exist at the given thermodynamic conditions. Will not be considered in the CPL solutions.')  
                        print('########################################################################################################################')  
                        #print('\n')  
                    else:                    
                        print('The following phase will co-exist everywhere within the CPL. Will not be considered in the CPL solutions.')  
                        print(self.inputs['sPhase'][s], sData['formEn'][s])     
                #at least one chem. spices is unconstrained 
                elif sceff_tot == 1  and (sData['formEn'][s] / sData['coeff'][s][0]) < pData['formEn'][p]: 
                    print('NOTE:')  
                    print('########################################################################################################################')  
                    print('Very poor (sumO_del_mu < del_Hf_primary) conditions for the following phase.')  
                    #print('The following phase will not be in the CPL solutions as their FE are poorer compared to FE_primary')  
                    print( self.inputs['pPhase'][p],": ",'Sec. phase: ', self.inputs['sPhase'][s],"; ", 'FE_sec (eV): ', round(sData['formEn'][s], 3),\
                             '; FE_prim (eV): ', round(pData['formEn'][p], 3)) #sData['var'][s])     
                    print('The above phase will not co-exist at the given thermodynamic conditions. Will not be considered in the CPL solutions.')  
                    print('########################################################################################################################')  
                #at least two chem. spices are unconstrained 
                elif sceff_tot == 2 and ((sData['formEn'][s] / sData['coeff'][s][0]) < pData['formEn'][p] or (sData['formEn'][s] / sData['coeff'][s][1]) < pData['formEn'][p]):
                    print('NOTE:')  
                    print('########################################################################################################################')  
                    print('Very poor (sumO_del_mu < del_Hf_primary) conditions for the following phase.')  
                    #print('The following phase will not be in the CPL solutions as their FE are poorer than FE_primary')  
                    print( self.inputs['pPhase'][p],": ",'Sec. phase: ', self.inputs['sPhase'][s],"; ", 'FE_sec (eV): ', round(sData['formEn'][s], 3),\
                             '; FE_prim (eV): ', round(pData['formEn'][p], 3)) #sData['var'][s])     
                    print('The above phase will not co-exist at the given thermodynamic conditions. Will not be considered in the CPL solutions.')  
                    print('########################################################################################################################')  
                #avoid non-equilibrium secondary phases; FE => 0
                elif sData['formEn'][s] >= 0:      
                    print('NOTE:')  
                    print('########################################################################################################################')  
                    print('Non-equilibrium conditions (sumOf_del_mu > 0) for the following phase.')  
                    print( self.inputs['pPhase'][p],": ",'Sec. phase: ', self.inputs['sPhase'][s],"; ", 'FE_sec (eV): ', round(sData['formEn'][s], 3),\
                             '; FE_prim (eV): ', round(pData['formEn'][p], 3)) #sData['var'][s])     
                    print('This above phase will not co-exist at the given thermodynamic conditions. Will not be considered in the CPL solutions.')  
                    print('########################################################################################################################')  
                #add sec phases FE equations to be solved with the primary counterpart  
                else:
                    sec.append(self.inputs['sPhase'][s])  #stoichiometry without constrain
                    sfe.append(sData['formEn'][s])        #constrained FE
                    #find the indices of chem. species of a sec phase in the primary-phase species list
                    #with Cu and S constrained, Cu2SnZnS2Se2 <==> SnZnSe2, in which Zn and Se have positions at 1 and 2
                    #Zn1Se1 matrix would be: [0, Zn, Se] = [0, 1, 1]; e.g., Zn1Se1: 1,1---->Zn,Se coefficients in Zn1Se
                    for i in range(sceff_tot):
                        svar  = sData['var'][s][i]
                        pIdx  = pData['var'][p].index(svar)
                        sIdx  = sData['var'][s].index(svar)
                        sceff = sData['coeff'][s][sIdx]
                        scfs[pIdx]   =  sceff
                        satms[pIdx]  =  svar
                    ssp.append(satms)
                    scoeff.append(scfs)
            self.leData['sPhase'].append(sec)
            self.leData['sSP'].append(ssp)
            self.leData['sceff'].append(scoeff)
            self.leData['sFE'].append(sfe)
                    
        return(self.leData)

    def solve_with_axes(self, leData, axisData, xaxis, yaxis):
        self.axSolData = collections.defaultdict(list)
        for p in range(len(leData['pPhase'])):
            lines, sp, pts, ax      = ([], [], [], [])
            self.axSolData['pPhase'].append(leData['pPhase'][p])
            self.axSolData['axis'].append(ax)
            self.axSolData['axis'][p].append(xaxis)
            self.axSolData['axis'][p].append(yaxis)
            pPhase = leData['pPhase'][p]
            diag   = axisData['axis'][p][2]
            didx   = leData['pSP'][p].index(diag)
            dceff  = leData['pceff'][p][didx]
            for s in range(len(leData['sPhase'][p])):
                sPhase = leData['sPhase'][p][s]
                em     = []
                for atm in range(len(leData['sSP'][p][s])):
                    atom  = leData['sSP'][p][s][atm]
                    #solve with the x-axis
                    #print(leData['pceff'][p])
                    #print(leData['sceff'][p][s])
                    if str(atom).isalpha() and atom == xaxis:
                        lines.append(leData['sPhase'][p][s])
                        sp.append(leData['sSP'][p][s])
                        #compare pPhase and sPhase slopes to see if the lines are parallel with yaxis = 0
                        mp  = leData['pceff'][p][2] / leData['pceff'][p][0]
                        ms  = leData['sceff'][p][s][2] / leData['sceff'][p][s][0]
                        #print('solve_axes: ', atom, leData['sPhase'][p][s])
                        if mp == ms:
                            #print('No solution with ', xaxis, '= 0: ', pPhase, 'and ', sPhase, ' are parallel')
                            #print('Solving with ', yaxis, ' = 0')
                            xA  = [leData['pceff'][p], leData['sceff'][p][s], axisData['axceff'][p][1][0]]
                            xB  = [leData['pFE'][p], leData['sFE'][p][s], axisData['axceff'][p][1][1]]
                            xA  = np.array(xA)
                            xB  = np.array(xB)
                            SOL, resid, rank, sigma  = linalg.lstsq(xA, xB)
                            #print(SOL)
                            SOL = SOL.tolist()
                            pts.append(SOL)
                        else:
                            #print('solving with ', yaxis, ' = 0')
                            xA  = [leData['pceff'][p], leData['sceff'][p][s], axisData['axceff'][p][0][0]]
                            xB  = [leData['pFE'][p], leData['sFE'][p][s], axisData['axceff'][p][0][1]]
                            xA  = np.array(xA)
                            xB  = np.array(xB)
                            #print(xA)
                            #print(xB)
                            SOL, resid, rank, sigma  = linalg.lstsq(xA, xB)
                            #when one of the obtained intersecting points lie outside CPL, it needs to be exptrapolated to x-/y-axis or cut to those axes
                            #check the indicative diagonal values: > 0 and < FE_constrained / diag_coefficient     
                            if SOL[2] > 0 or SOL[2] < leData['pFE'][p]/dceff:
                                #print('one intersecting point is outside the CPL with ', yaxis, ' = 0')
                                #print('resetting the solving axis,', xaxis, ' = 0')
                                xAfit  = [leData['pceff'][p], leData['sceff'][p][s], axisData['axceff'][p][1][0]]
                                xBfit  = [leData['pFE'][p], leData['sFE'][p][s], axisData['axceff'][p][1][1]]
                                xAfit  = np.array(xAfit)
                                xBfit  = np.array(xBfit)
                                SOLfit, resid, rank, sigma  = linalg.lstsq(xAfit, xBfit)
                                #print(SOLfit)
                                SOL = SOLfit.tolist()
                                pts.append(SOL)
                            else:
                                #print(SOL)
                                SOL = SOL.tolist()
                                pts.append(SOL)
                    if str(atom).isalpha() and atom == yaxis:
                        lines.append(leData['sPhase'][p][s])
                        sp.append(leData['sSP'][p][s])
                        #compare pPhase and sPhase slopes to see if the lines are parallel with xaxis = 0
                        mp  = leData['pceff'][p][2] / leData['pceff'][p][1]
                        ms  = leData['sceff'][p][s][2] / leData['sceff'][p][s][1]
                        #print('solve_axes: ', atom, leData['sPhase'][p][s])
                        if mp == ms:
                            #print('No solution with ', xaxis, '= 0: ', pPhase, 'and ', sPhase, ' are parallel')
                            #print('Solving with ', yaxis, ' = 0')
                            xA  = [leData['pceff'][p], leData['sceff'][p][s], axisData['axceff'][p][0][0]]
                            xB  = [leData['pFE'][p], leData['sFE'][p][s], axisData['axceff'][p][0][1]]
                            xA  = np.array(xA)
                            xB  = np.array(xB)
                            SOL, resid, rank, sigma  = linalg.lstsq(xA, xB)
                            #print(SOL)
                            SOL = SOL.tolist()
                            pts.append(SOL)
                        else:
                            #print('solving with ', xaxis, ' = 0')
                            xA  = [leData['pceff'][p], leData['sceff'][p][s], axisData['axceff'][p][1][0]]
                            xB  = [leData['pFE'][p], leData['sFE'][p][s], axisData['axceff'][p][1][1]]
                            xA  = np.array(xA)
                            xB  = np.array(xB)
                            #print(xA)
                            #print(xB)
                            SOL, resid, rank, sigma  = linalg.lstsq(xA, xB)
                            #when one of the obtained intersecting points lie outside CPL, it needs to be exptrapolated to x-/y-axis or cut to those axes
                            #check the indicative diagonal values: > 0 and < FE_constrained / diag_coefficient     
                            if SOL[2] > 0 or SOL[2] < leData['pFE'][p]/dceff:
                                #print('one intersecting point is outside the CPL with ', xaxis, ' = 0')
                                #print('resetting the solving axis,', yaxis, ' = 0')
                                xAfit  = [leData['pceff'][p], leData['sceff'][p][s], axisData['axceff'][p][0][0]]
                                xBfit  = [leData['pFE'][p], leData['sFE'][p][s], axisData['axceff'][p][0][1]]
                                xAfit  = np.array(xAfit)
                                xBfit  = np.array(xBfit)
                                SOLfit, resid, rank, sigma  = linalg.lstsq(xAfit, xBfit)
                                #print(SOLfit)
                                #pts[s].append(SOLfit.tolist())
                                SOL = SOLfit.tolist()
                                pts.append(SOL)
                            else:
                                #print(SOL)
                                #pts[s].append(SOL.tolist())
                                SOL = SOL.tolist()
                                pts.append(SOL)
            self.axSolData['points'].append(pts)
            self.axSolData['sPhase'].append(lines)
            self.axSolData['sp'].append(sp)
        return(self.axSolData)

    def solve_with_diag(self, leData, axisData, solAxis):
        self.diaSolData = collections.defaultdict(list)
        for p in range(len(leData['pPhase'])):
            xaxis, yaxis, diag  = (axisData['axis'][p][0], axisData['axis'][p][1], axisData['axis'][p][2])
            lines, sp, pts, ax      = ([], [], [], [])
            self.diaSolData['pPhase'].append(leData['pPhase'][p])
            self.diaSolData['axis'].append(ax)
            self.diaSolData['axis'][p].append(diag)
            for s in range(len(leData['sPhase'][p])):
                em   = []
                for atm in range(len(leData['sSP'][p][s])):
                    atom  = leData['sSP'][p][s][atm]
                    #solve with the diagonal
                    if str(atom).isalpha() and (atom == xaxis or atom == yaxis):
                        lines.append(leData['sPhase'][p][s])
                        #print('solve_axes: ', solAxis, leData['sPhase'][p][s])
                        diaA = [leData['pceff'][p],  leData['sceff'][p][s], axisData['axceff'][p][2][0]]
                        diaB = [leData['pFE'][p], leData['sFE'][p][s], axisData['axceff'][p][2][1]]
                        diaA = np.array(diaA)
                        diaB = np.array(diaB)
                        #print(diaA)
                        #print(diaB)
                        SOLd, resid, rank, sigm  = linalg.lstsq(diaA, diaB)
                        #print(SOLd)
                        pts.append(SOLd.tolist())
            self.diaSolData['sPhase'].append(lines)
            self.diaSolData['sp'].append(sp)
            self.diaSolData['points'].append(pts)
        return(self.diaSolData)

    def write_cpl_data(self, cplDataXY, cplDataD):
        with open('cpl.dat', 'w') as f:
            f.write(str('#') + str('chemical potentials for the given thermodynamic conditions'))
            f.write("\n")
            for p  in range(len(cplDataXY['pPhase'])):
                xaxis  = cplDataXY['axis'][p][0]
                yaxis  = cplDataXY['axis'][p][1]
                diag   = cplDataD['axis'][p][0]
                pPhase = cplDataXY['pPhase'][p]
                sCounted = []
                #print(xaxis, yaxis, diag)
                f.write(str('#') + str(pPhase))
                f.write("\n")
                f.write(str('#') + str('axis: ') + "\t" + str(xaxis) + "\t" +str(yaxis)+ "\t" + str(diag))
                f.write("\n")
                for s  in range(len(cplDataXY['sPhase'][p])):
                    if s < (len(cplDataXY['sPhase'][p]) - 1) and cplDataXY['sPhase'][p][s] == cplDataXY['sPhase'][p][s+1]:
                        sPhaseXY = cplDataXY['sPhase'][p][s]
                        sCounted.append(sPhaseXY) 
                        xval1  = round(cplDataXY['points'][p][s][0], 3)
                        xval2  = round(cplDataXY['points'][p][s+1][0], 3)
                        yval1  = round(cplDataXY['points'][p][s][1], 3)
                        yval2  = round(cplDataXY['points'][p][s+1][1], 3)
                        zval1  = round(cplDataXY['points'][p][s][2], 3)
                        zval2  = round(cplDataXY['points'][p][s+1][2], 3)
                        f.write(str(sPhaseXY) + "\t" + str(xval1) + "\t" +str(yval1)+ "\t" + str(zval1))
                        f.write("\n")
                        f.write(str(sPhaseXY) + "\t" + str(xval2) + "\t" +str(yval2)+ "\t" + str(zval2))
                        f.write("\n")
                    else: 
                        sPhase = cplDataXY['sPhase'][p][s]
                        if sPhase not in sCounted: 
                            xval1  = round(cplDataXY['points'][p][s][0], 3)
                            xval2  = round(cplDataD['points'][p][s][0], 3)
                            yval1  = round(cplDataXY['points'][p][s][1], 3)
                            yval2  = round(cplDataD['points'][p][s][1], 3)
                            zval1  = round(cplDataXY['points'][p][s][2], 3)
                            zval2  = round(cplDataD['points'][p][s][2], 3)
                            f.write(str(sPhase) + "\t" + str(xval1) + "\t" +str(yval1)+ "\t" + str(zval1))
                            f.write("\n")
                            f.write(str(sPhase) + "\t" + str(xval2) + "\t" +str(yval2)+ "\t" + str(zval2))
                            f.write("\n")

