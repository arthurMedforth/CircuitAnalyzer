# -*- coding: utf-8 -*-
"""
@author: am3235
@Candidate number: 11071
@Date: 20/04/2020
This is a program to read-in, analyse and write-out electrical circuit information described in an input file.
"""
### Import libraries
import numpy as np
import sys
import os
import cmath
import math


###########################################################################################################################################################################################
######################################## FUNCTIONS ########################################################################################################################################
###########################################################################################################################################################################################


########## Functions to multiply matrices within list #################################################################
def matMult_1(List_of_matrices):
    '''
    :param List_of_matrices: List of lists (matrices) to be multiplied together

    :return: List of matrices: Return the single element list of lists

    Operation:
    Take input argument and multiply all lists inside in order
    '''

    if len(List_of_matrices)==1: # If there is only one element, no multiplication can be done
        return List_of_matrices  # so return it
    else:
        Result = np.dot(List_of_matrices[0],List_of_matrices[1])  # Matrix dot product of the first two elements
        List_of_matrices.pop(1)  # Delete the second element (index = 1)
        List_of_matrices.pop(0)  # and then the first (index = 0)
        List_of_matrices.insert(0,Result) # Insert result of multiplcation in place of two deleted matrices
        return matMult_1(List_of_matrices) # Recall the same function

def matMult(matrixList,listLength):
    '''
    :param matrixList: List of matrices
    :param listLength: Length of list

    :return: matrixList1: Replication of matrixList

    Operation:
    Multiply all matrices within matrixList by eachother in order
    '''

    matrixList1=matrixList[:]      # Replicate to avoid global errors
    for i in range(listLength):
     if len(matrixList1) == 1:
        break
     elif len(matrixList1) >= 2:
        f = np.dot(matrixList1[0], matrixList1[1])
        matrixList1.pop(1)
        matrixList1.pop(0)
        matrixList1.insert(0,f)
        matMult(matrixList1,listLength)
    return matrixList1

############### Function to find type of circuit #######################################################################
def whatCircuit(LVector,CVector,resVector):
    '''
    :param LVector: List of any inductor values in the circuit
    :param CVector: List of any capacitor values in the circuit
    :param resVector: List of any resistor values in the circuit

    :return: R, LR, RC,LC, LCR
    Corresponding to different types of circuit

    Operation:
    Return type of circuit (R,LR,RC,LC or LCR)
    '''
    ### Initialise boolean variables
    CinCirc=bool
    LinCirc=bool
    RinCirc=bool
    LCR=bool
    RC=bool
    LR=bool
    LC=bool
    RCirc=bool

    # Loop through each list individually
    for o in range(len(node1Vector)):
        if CVector[o] != 0:
            CinCirc = True  # Set CinCirc bool if capacitors present
            break

    for o in range(len(node1Vector)):
        if LVector[o] != 0:
            LinCirc = True  # Set LinCirc bool if capacitors present
            break

    for o in range(len(node1Vector)):
        if resVector[o] != 0:
            RinCirc = True # Set RinCirc bool if capacitors present
            break

    if LinCirc == True and CinCirc == True and RinCirc == True: # Circuit contains all L, R and C
        LCR = True                                              # LCR type
    elif LinCirc == True and CinCirc == True: # Circuit contains all and C
        LC = True                             # LC type
    elif LinCirc == True:          # Circuit contains all L and R
        LR = True                  # LR type
    elif CinCirc == True:   # Circuit contains all C and R
        RC = True           # RC type
    else:
        RCirc = True # Circuit just contains R so RCirc type

    return LCR, LR, LC, RC, RCirc

############ Function to find indices of inductive elements ###########################################################
def getIndexPositionsL(listOfElements):
 '''
 :param listOfElements: branches vector for whole circuit

 :return: indVecL: list of indices of inductive branches

 Operation:
 Return indices of inductive branches within branches
 '''

 indVecL = [] # Initialise output list
 for i in range(len(listOfElements)): # Loop through branches
     if listOfElements[i][4] == 2: # Look for '2' tag indicating inductor
         indVecL.append(i) # If found - append the index into indVecL
     if i == len(listOfElements)-1: # If checked final branch in branches
        return(indVecL) # return the list of indices
     else:
         pass

############ Function to find indices of capacitive elements ###########################################################
def getIndexPositionsC(listOfElements):
 '''
 :param listOfElements: branches vector for whole circuit

 :return: indVecL: list of indices of inductive branches

 Operation:
 Return indices of capacitive branches within branches
 '''

 indVecC = [] # Initialise output list
 for i in range(len(listOfElements)):
     if listOfElements[i][3]==1: # Look for '1' tag indicating capacitor
         indVecC.append(i) # If found - append the index into indVecL
     if i == len(listOfElements)-1: # If checked final branch in branches
        return(indVecC) # return the list of indices
     else:
         pass

############ Function to find ABCD matrices of RC circuits #############################################################
def capFunc(ABCDmatrices, indexPosList, capBranch, frequency,fEnd):
    '''
    :param ABCDmatrices: Holds ABCD matrices for all the resistive components
    :param indexPosList: List of capacitive branch indices
    :param capBranch: All capacitive branches from branches list
    :param frequency: List of frequencies to be used in impedance calculation
    :param fEnd: Final frequency from frequency vector

    :return: CascTFmatrices: contatins all cascade ABCD matrices for the different frequencies

    Operation:
    For each frequency, capFunc creates and inserts the series or shunt capacitive ABCD matrices into their appropriate
    positions using the index vector. At the end of the each frequency, these same elements are deleted from the
    ABCDmatrices list so that the process can be repeated with a new frequency. It outputs cascTFmatrices which holds
    all n-frequency cascade ABCD matrices for the circuit.
    '''
    ### Initialise lists
    capMat = []
    cascTF = []
    cascTFmatrices = []

    for l in frequency: # Loop through all frequencies
        for m in range(len(capBranch)): # Loop through each of the capacitor branches in capBranch list
            indNeed=indexPosList[m]  # Index for inserting the capacitor ABCD matrix
            if capBranch[m][1] == 0: # Shunt impedance due to node 2 = 0
                capMat = [[1, 0], [2 * np.pi * l * capBranch[m][2]*1j, 1]] # Create ABCD matrix
                ABCDmatrices.insert(indNeed, capMat)  # Insert matrix into ABCD matrices at appropriate position.
                if m == len(capBranch)-1: # If this is last capacitor branch
                    cascTF = matMult(ABCDmatrices, len(ABCDmatrices)) # Do the matrix multiplication of the ABCDmatrices
                    cascTFmatrices.append(cascTF) # Append result into cascTFmatrices
                    for index in sorted(indexPosList, reverse=True):  # delete capacitive elements from ABCDmatrices
                        del ABCDmatrices[index]  # must do in reverse order so that indexing works
                else:  # Otherwise continue onto next capBranch element
                    pass
            else: # Series impedance
                capMat = [[1, (1 / (2 * np.pi * l * capBranch[m][2]*1j))], [0, 1]]  # Create ABCD matrix
                ABCDmatrices.insert(indNeed, capMat) # Insert it into ABCDmatrices at appropriate position
                if m == len(capBranch)-1: # If the last branch in capBranch
                    cascTF = matMult(ABCDmatrices, len(ABCDmatrices)) # Execute multiplication
                    cascTFmatrices.append(cascTF) # Append result into cascTFmatrices
                    for index in sorted(indexPosList, reverse=True):  # delete capacitive elements from ABCDmatrices
                        del ABCDmatrices[index]  # must do in reverse order so that indexing works
                else:  # Otherwise continue onto next capBranch element
                    pass
    return(cascTFmatrices)

############ Function to find ABCD matrices of LR circuits #############################################################
def LFunc(ABCDmatrices, indexPosListL, LBranch, frequency,fEnd):
    '''
    :param ABCDmatrices: Holds ABCD matrices for all the resistive components
    :param indexPosListL: List of inductive branch indices
    :param LBranch: All inductive branches from branches list
    :param frequency: List of frequencies to be used in impedance calculation
    :param fEnd: Final frequency from frequency vector

    :return: CascTFmatrices: contatins all cascade ABCD matrices for the different frequencies

    Operation:
    For each frequency, LFunc creates and inserts the series or shunt inductive ABCD matrices into their appropriate
    positions using the index vector. At the end of the each frequency, these same elements are deleted from the
    ABCDmatrices list so that the process can be repeated with a new frequency. It outputs cascTFmatrices which holds
    all n-frequency cascade ABCD matrices for the circuit.
    '''
    ### Initialise lists
    LMat = []
    cascTF = []
    cascTFmatrices = []
    lenindexPosList=len(indexPosListL)

    for l in frequency: # Loop through frequencies
        for m in range(len(LBranch)): # Loop though each branch in LBranch
            indNeed=indexPosListL[m] # Index needed for insertion
            if LBranch[m][1] == 0: # Shunt impedance due to node 2 = 0
                LMat = [[1, 0], [(1 / (2 * np.pi * l * LBranch[m][2])*1j), 1]] # Create ABCD matrix
                ABCDmatrices.insert(indNeed, LMat) # Insert into ABCD matrices at appropriate index
                if m == len(LBranch)-1:  # If last branch used
                    cascTF = matMult(ABCDmatrices, len(ABCDmatrices)) # Execute multiplication of ABCDmatrices
                    cascTFmatrices.append(cascTF) # Append cascade ABCD matrix into cascTFmatrices
                    for index in sorted(indexPosListL, reverse=True):  # delete capacitive elements from branches
                        del ABCDmatrices[index]  # must do in reverse order so that indexing works
                else: # Otherwise continue onto next branch
                    pass
            else: # Series impedance
                LMat = [[1, (2 * np.pi * l * LBranch[m][2])*1j], [0, 1]] # Create ABCD matrix
                ABCDmatrices.insert(indNeed,LMat) # Insert inductive ABCD matrix at appropriate index
                if m == len(LBranch)-1:  # If last branch used
                    cascTF = matMult(ABCDmatrices, len(ABCDmatrices)) # Execute multiplication of ABCDmatrices
                    cascTFmatrices.append(cascTF) # Append cascade ABCD matrix into cascTFmatrices
                    for index in sorted(indexPosListL, reverse=True):  # delete capacitive elements from branches
                        del ABCDmatrices[index]  # must do in reverse order so that indexing works
                else: # Otherwise continue onto next branch
                    pass
    return(cascTFmatrices)

############ Function to find ABCD matrices of LCR or LC circuits ######################################################
def LCRFunc(ABCDmatrices,ABCDmatricesInit, indexPosListL,indexPosListC,LBranch,capBranch, frequency,fEnd):
    '''
    :param ABCDmatrices: Holds ABCD matrices for all the resistive components
    :param ABCDmatricesInit: Replica of ABCD matrices
    :param indexPosListL: Index position list for inductive elements
    :param indexPosListC: Index position list for capacitive elements
    :param LBranch: inductive branches list
    :param capBranch: capacitive branches list
    :param frequency: List of frequencies to be used in impedance calculation
    :param fEnd: Final frequency from frequency vector

    :return: cascTFmatrices: List containing all the cascade ABCD matrices

    Operation:
    For each frequency, LCRFunc creates and inserts the series or shunt, inductive or capacitive, ABCD matrices into
    their appropriate positions using the index vectors. At the end of the each frequency, delete these same elements
    from the ABCDmatrices list so that the process can be repeated with a new frequency. It outputs cascTFmatrices
    which holds all n-frequency cascade ABCD matrices for the circuit.
    '''
    ### Initialise lists and variables
    LMat = []
    capMat = []
    cascTF = []
    cascTFmatrices = []
    indexDeleteVect= indexPosListC + indexPosListL # Concatenate index vectors
    indexDeleteVect.sort() # and sort into ascending order

    for l in frequency:
        ABCDmatrices=ABCDmatricesInit[:] # Create replica ABCDmatrices list for multiplication
        for m in range(len(LBranch)): # Loop though inductive branches
            indNeedL=indexPosListL[m] # Index for insertion into ABCDmatrices
            if LBranch[m][1] == 0: # Shunt component
                LMat = [[1, 0], [(1 / (2 * np.pi * l * LBranch[m][2]*1j)), 1]] # Create ABCD matrix
                ABCDmatrices.insert(indNeedL, LMat) # Insert into ABCD matrices
            else: # Series component
                LMat = [[1, (2 * np.pi * l * LBranch[m][2])*1j], [0, 1]]
                ABCDmatrices.insert(indNeedL,LMat) # Insert into ABCD matrices

        for m in range(len(capBranch)): #Loop through capacitive branches
            indNeedC = indexPosListC[m] # Index for insertion into ABCDmatrices
            if capBranch[m][1] == 0: # Shunt component
                capMat = [[1, 0], [(2 * np.pi * l * capBranch[m][2]*1j), 1]] # Create ABCD matrix
                ABCDmatrices.insert(indNeedC, capMat) # Insert into ABCDmatrices
                if m == len(capBranch)-1: # If last capacitive branch
                    ABCDmodel=ABCDmatrices[:] # Create replica of ABCD matrices
                    cascTF = matMult_1(ABCDmodel) # Execute multiplication
                    cascTFmatrices.append(cascTF) # Append cascade ABCD matrix into ABCDmatrices
                else:
                    pass
            else: # Series component
                capMat = [[1, (1 / (2 * np.pi * l * capBranch[m][2]*1j))], [0, 1]] # Create ABCD matrix
                ABCDmatrices.insert(indNeedC, capMat) # Insert ABCDmatrix at appropriate index
                if m == len(capBranch)-1: # If the final branch
                    ABCDmodel=ABCDmatrices[:] # Create replica
                    cascTF = matMult_1(ABCDmodel) # Execute multiplication
                    cascTFmatrices.append(cascTF) # Append result into cascTFmatrices
                else: # Otherwise, continue onto next branch
                    pass

    return(cascTFmatrices)


######### Function to find delimiters ##################################################################################
def find_string(string_to_search, name):
    '''
    :param string_to_search: the text string to be searched
    :param name: text string of the item to be sought in string_to_search

    :return: Found: 1 if found, 0 if not

    Operation:
    find_string looks through lines of the input file to find a string given by 'name' parameter
    '''

    pp=string_to_search.find(name) # Look for string with in-house 'find'
    if (pp > -1):  # name has been found
        Found = 1
        return(Found)
    else:
        Found = 0
        return(Found)

########## Function to find values in circuit branches and terminal information ########################################
def find_float(string_to_search,name,noisy):
    '''
    :param string_to_search: the text string to be searched
    :param name:    text string of the item to be sought in string_to_search
    :param noisy:  boolean, if true produces more output to the screen

    :return: rtn - integer value found, or particular negative values when not found
              ok - boolean, true if name and integer found otherwise false

    Operation:
    find_float searches for a string given by name in the input file line, looks for an equals sign and attempts to
    return the float type number after the sign.
    '''

    pp=string_to_search.find(name)    # looks for name in string_to_search
    if (pp>-1): # name has been found
        lname=len(name)
        idx=pp+lname   # points at character after name
        new_str=string_to_search[idx:]  # copy from idx to end
        qqe=new_str.find('=')  # search for equals sign
        if (qqe>-1):   # equals found
            idx=qqe+1  # point at character after equals sign
            valstring=new_str[idx:]
            substring=valstring.split(' ') #variables need to be separated by spaces
            teststring=substring[0].strip() # strip gets rid of spare spaces
            try:
                rtn=float(teststring)  #convert string to integer
                ok=True
            except ValueError:  # Right type invalid value
                typeTString=type(teststring)
                if typeTString is str: # If the supposed number is a string
                    try: # Open file
                        W = open(sys.argv[2], "w+") # If called from the command line
                    except:
                        W = open('My_Circuit.dat', "w+") # Otherwise open file specified in run parameters ^^^
                    W.close() # Close the blank file
                    sys.exit('Input formatting error') # Leave with error message
                else:
                    rtn=0 # Return zero otherwise
                    ok=False
        else:  # equals not found
            ok=False
            rtn=0
    else:  # name not found in string_to_search
        rtn=0
        ok = False
    return rtn, ok

##############################################################################################################################################################################################################################
###########################################  MAIN PROGRAM ####################################################################################################################################################################
##############################################################################################################################################################################################################################

# Initialise variables and lists
node1Vector=[]
node2Vector=[]
resVector=[]
GVector=[]
branches=[]
ABCDmatrices=[]
ABCDmatrix=[]
VTVec=[]
termsVector = []
outputVector = []
freqVect=[]
CVector=[]
LVector=[]
cCheck=[]
lCheck=[]
Lfrequency=[]
capBranch=[]
LBranch=[]
outqs=[]
ZinV=[]
VinV=[]
IinV=[]
PinV=[]
inputMatV=[]
outputMV=[]
VoutV=[]
IoutV=[]
PoutV=[]
ZoutV=[]
AvV=[]
AiV=[]
AidBV=[]
outV = []
LFend = 0
LFstart = 0
indexPos = 0
LFREQ = bool
FREQ = bool


if sys.argv[1] == '-i': # Check if Time response program called from command line
    input_filename = sys.argv[2] + '.net' # if so, take the third argument and concatentate .net onto end
    pointsIFFT = 256 # Set bit precision of IFFT: 2^n
else:
    input_filename = sys.argv[1] # Otherwise, take second argument as the input filename

### Display valuable program input parameters
print("This is the name of the script: ", sys.argv[0])
print("Number of arguments: ", len(sys.argv))
print("The arguments are: " , str(sys.argv))
if (len(sys.argv)<2): # If not enough arguments produce error
    print("\n\tError, program needs two arguments to run\n" )
    sys.exit(1) # Exit program

# Open file
try:
    fin=open(input_filename,'rt') # File specified in command line or other
except:
    print('File <%s> not found'%(input_filename))
    current_location=os.getcwd()
    # Prints working directory
    print("executing program in directory: "+current_location)
    sys.exit('Check spelling and/or location of file')  # exits the program, returning a value of 1 to the operating system

# Now read file
file_lines=fin.readlines()

for index,line in enumerate(file_lines): # Loop through lines of the input file
    if line.find('#') > -1:  # If there is a # at start of line, the information is a comment
        continue  # Skip to next line
    else:
        ### Look for all parameters with numeric value on each line
        ### If possible, retrieve float value and output as much information as possible
        n1, n1_found=find_float(line,"n1",True)
        n2, n2_found=find_float(line,"n2",True)
        R, res_found=find_float(line,"R",True)
        G, ad_found=find_float(line,"G",True)
        VT, vSource_found=find_float(line,"VT",True)
        IN, iSource_found=find_float(line,"IN",True)
        RS, RS_found=find_float(line,"RS",True)
        GS, GS_found=find_float(line,"GS",True)
        RL, RL_found=find_float(line,"RL",True)
        C, C_found=find_float(line,"C",True)
        L, L_found=find_float(line,"L",True)
        fStart, fStart_found=find_float(line,'Fstart',True)
        fEnd, fEnd_found=find_float(line,'Fend',True)
        Nfreqs, Nfreqs_found=find_float(line,'Nfreqs',True)
        LFstart, LFstart_found = find_float(line, 'LFstart', True)
        LFend, LFend_found = find_float(line, 'LFend', True)

        ### Look for OUTPUT delimiters
        delimiterstartFound = find_string(line, "<OUTPUT>")
        delimiterendFound = find_string(line,"</OUTPUT>")


    if delimiterstartFound==1: # Delimiter found
        startIndex=index+1   # Save start index
    elif delimiterendFound==1: # Delimiter found
        endIndex=index+1  # Save end index

    if LFend != 0 and LFstart != 0: # If logarithmic frequency information found
        LFREQ = True # Set boolean variables
        FREQ = False
        LFend = int(LFend)       #Convert to integers for input into numpy linspace
        LFstart = int(LFstart)
        Nfreqs = int(Nfreqs)
        Lfrequency = np.geomspace(LFstart, LFend, Nfreqs) # Produce frequency vector using numpy.linspace

    elif fEnd !=0 and fStart!=0 and Nfreqs!=0 and LFend == 0: # If input file asks for normal frequency sweep
       FREQ = True   # Set boolean variables
       LFREQ = False
       frequency = [fStart]
       Nfreqs = int(Nfreqs)
       fGap = (fEnd-fStart)/(Nfreqs-1)

       if sys.argv[1] == '-i':  # If time response called from command line
            freqRangeIFFT = np.linspace(-fEnd/2,fEnd/2,pointsIFFT) # produce IFFT frequency range with frequency information and IFFT points
       else:
            pass

       for f in range(1,Nfreqs):     # Loop through and calculate frequency vector
           frequency.append(fStart+(fGap*f))

    else: # if no frequency information found, continue to next line
        pass

    if VT==0 and IN==0:  # If no source information found, move to next line
        pass
    else:
        termsVector.append([VT,IN,RS,GS]) # Append any found terms values found into termsVector
        ### Convert to Thevenin equivalent circuit for uniform procedure in subsequent section
        if termsVector[0][0] == 0 and termsVector[0][1] != 0 and termsVector[0][2] != 0: # If Norton current given with RS...
            termsVector[0][0] = termsVector[0][1] * termsVector[0][2] # calculate Thevenin source
        elif termsVector[0][0] == 0 and termsVector[0][1] != 0 and termsVector[0][3] != 0: # If Norton current given with GS...
            termsVector[0][0] = termsVector[0][1] * 1/termsVector[0][3] # Calculate VT
            termsVector[0][2] = 1/termsVector[0][3]  # Calculate RS

    if RL!=0: # If found...
        termsVector[0].append(RL) # Append RL into termsVector

    if n1 == 0 and n2 ==0: # If nodes 1 and 2, pass to next line
        pass
    else: # If not, append their values into lists
        n1Int=int(n1)
        node1Vector.append(n1Int)
        n2Int=int(n2)
        node2Vector.append(n2Int)
        resVector.append(R) # Also append any resistance
        GVector.append(G)   # or admittance
        if C!=0: # If capacitance value not equal to zero...
            CVector.append(C) # Append capacitor values to capacitor vector
            LVector.append(0) # Append zero to inductor vector
        elif L!=0: # If inductance value not equal to zero
            LVector.append(L) # Append inductor value into inductor vector
            CVector.append(0) # Append zero into capacitor vector
        else:
            CVector.append(0) # Otherwise append zeros into both
            LVector.append(0)

for i in range(len(node1Vector)):
    if node1Vector[i] == 0: # Node 1 value could be zero and node 2 non-zero
        node1Vector[i] = node2Vector[i] # if so...
        node2Vector[i] = 0 # swap node 1 and node 2 around

if termsVector[0][0] != 0 and termsVector[0][3] != 0: # Check again to see if termsVector is Thevening equivalent circuit
    termsVector[0][2] = 1/termsVector[0][3] # Change admittance to resistance if needed
else:
    pass

for i in range(startIndex, endIndex-1):
    lineOI=file_lines[i]
    commentFound = find_string(lineOI,"#") # Look for # to determine if a comment line
    if commentFound == 0: # Comment line found
        outputVector.append(lineOI)                         # Append output vector with output quantities needed
        for j in range(len(outputVector)):                  # Remove '\n' produced by new line in file script interpretation
            outputVector[j] = outputVector[j].rstrip("\n")  # Strip output quantity strings of any python new line or
            outputVector[j] = outputVector[j].rstrip("\r")  # carriage return characters
    else:
        pass

ABCDimpedanceVector=resVector  # Initialise an ABCDimpedance vector by equating it to resVector
for m in range(len(node1Vector)): # For loop with range of node1Vector
    cCheck.append(0) # Initialise cCheck to all zeros
    lCheck.append(0) # Initialise lCheck to all zeros

for i in range(len(node1Vector)): # For loop with range of node1Vector
    if node2Vector[i] == 0: # Shunt connection
        if resVector[i] == 0 and GVector[i] != 0: # only admittance given
            resVector[i]=1/GVector[i] # Get resistance from admittance
            cCheck[i]=0 # No capacitor or inductor...
            lCheck[i]=0 # so these can be set to zero
        elif CVector[i] != 0:
            ABCDimpedanceVector[i]=CVector[i] # Append capacitor value into ABCDimpedance vector
            cCheck[i]=1 # Set tags in cCheck and
            lCheck[i]=0 # lCheck
        elif LVector[i]!= 0:
            ABCDimpedanceVector[i]=LVector[i] # Append inductor value into ABCDimpedance vector
            lCheck[i]=2 # Set tags in lCheck and
            cCheck[i]=0 # cCheck
        else:
            pass
    else: # Series connection
        if resVector[i]==0 and GVector[i] != 0: # If admittance given...
            resVector[i]=1/GVector[i] # calculate resistance
            cCheck[i]=0 # Set tags in cCheck and
            lCheck[i]=0 # lCheck
        elif CVector[i] != 0:
            ABCDimpedanceVector[i] = CVector[i] # Append capacitor value into ABCDimpedanceVector
            cCheck[i] = 1 # Set tags in cCheck and
            lCheck[i] = 0 # lCheck
        elif LVector[i] != 0:
            ABCDimpedanceVector[i] = LVector[i] # Append inductor value into ABCDimpedanceVector
            lCheck[i] = 2 # Set tags in lCheck and
            cCheck[i] = 0 # cCheck

for k in range(len(node1Vector)): # Now go through ABCDimpedance vector
    if node2Vector[k] == 0: # If shunt connection
        if CVector[k] == 0 and LVector[k] == 0: # and capacitor and inductor values are also zero
            ABCDimpedanceVector[k]=1/resVector[k] # change ABCDimpedanceVector to admittance using the resistance
    else: # Otherwise leave ABCDimpedanceVector as is
        pass # Move onto to next node connection

num_branches=len(node1Vector)    # Number of node connections
for p in range(num_branches): # Loop through all node connections
        branches.append([node1Vector[p],node2Vector[p],ABCDimpedanceVector[p],cCheck[p],lCheck[p]]) # Create list of lists of node numbers and impedances
        if p == (num_branches-1):                                                                     # or capacitances or inductances for ABCD matrix calculation
            branches.sort() #sort branches list for

indexPosListC = getIndexPositionsC(branches) # Get indices of all capacitive branches
indexPosListL = getIndexPositionsL(branches) # Get indices of all inductive branches
indicesLCR = indexPosListC + indexPosListL # Concatenate the lists of indices
indicesLCR.sort() # Sort into ascending order

LCR,LR,LC,RC,RCirc = whatCircuit(LVector,CVector,resVector) # Return type of circuits

if RC == True:
    for item in indexPosListC: # for each relevant index
        capBranch.append(branches[item]) # Append any capacitive elements of branches into capBranch
    for index in sorted(indexPosListC, reverse = True): # delete capacitive elements from branches
        del branches[index]   # must do in reverse order so that indexing works

if LR == True:
    for item in indexPosListL: # for each relevant index
        LBranch.append(branches[item]) # Append any inductive elements of branches into LBranch
    for index in sorted(indexPosListL, reverse = True): # delete capacitive elements from branches
        del branches[index]   # must do in reverse order so that indexing works
elif LCR == True or LC == True:
    for item in indexPosListL: # for each relevant index
        LBranch.append(branches[item]) # Append any inductive elements of branches into LBranch
    for item in indexPosListC: # for each relevant index
        capBranch.append(branches[item]) # Append any capacitive elements of branches into capBranch
    for index in sorted(indicesLCR, reverse = True): # delete capacitive elements from branches
        del branches[index]   # must do in reverse order so that indexing works
else:
    pass

if branches == []:
    pass
else:
    for l in range(len(branches)): # For the resistive node connections or branches, create appropriate ABCDmatrix
        if branches[l][1]==0: # Shunt connection
            ABCDmatrix = [[1,0],[branches[l][2],1]]
            ABCDmatrices.append(ABCDmatrix) # Append matrix into ABCDmatrices
        else:  #Series connection
            ABCDmatrix = [[1,branches[l][2]],[0,1]]
            ABCDmatrices.append(ABCDmatrix) # Append matrix into ABCDmatrices


### Calculate cascTFmatrices vector for time response
if sys.argv[1] == '-i': # If time response called from the command line
    # Use IFFT frequency range in cascTFmatrices creation function
    if RC == True:
       cascTFmatrices=capFunc(ABCDmatrices, indexPosListC,capBranch,freqRangeIFFT,fEnd)
    elif LR == True:
       cascTFmatrices = LFunc(ABCDmatrices, indexPosListL, LBranch, freqRangeIFFT, fEnd)
    elif LCR == True or LC == True:
       ABCDinit = ABCDmatrices
       cascTFmatrices = LCRFunc(ABCDmatrices,ABCDinit,indexPosListL,indexPosListC,LBranch,capBranch,freqRangeIFFT,fEnd)
    elif RCirc == True:
        length4index = len(ABCDmatrices)
        cascTF = matMult_1(ABCDmatrices)

### Calculate cascTFmatrices vector for frequency or logarithmic frequency sweep
else:
    if RC == True:
        # Calculate cascTFmatrices using the required frequency sweep - determined by boolean variables
        if FREQ == True:
            cascTFmatrices=capFunc(ABCDmatrices, indexPosListC,capBranch,frequency,fEnd)
        elif LFREQ == True:
            cascTFmatrices = capFunc(ABCDmatrices, indexPosListC, capBranch, Lfrequency, fEnd)
    elif LR == True:
        # Calculate cascTFmatrices using the required frequency sweep - determined by boolean variables
        if FREQ == True:
           cascTFmatrices = LFunc(ABCDmatrices, indexPosListL, LBranch, frequency, fEnd)
        elif LFREQ == True:
           cascTFmatrices = LFunc(ABCDmatrices, indexPosListL, LBranch, Lfrequency, fEnd)
    elif LCR == True or LC == True:
        # Calculate cascTFmatrices using the required frequency sweep - determined by boolean variables
        ABCDinit = ABCDmatrices
        if FREQ == True:
           cascTFmatrices = LCRFunc(ABCDmatrices,ABCDinit,indexPosListL,indexPosListC,LBranch,capBranch,frequency,fEnd)
        elif LFREQ == True:
           cascTFmatrices = LCRFunc(ABCDmatrices, ABCDinit, indexPosListL, indexPosListC, LBranch, capBranch, Lfrequency, fEnd)
    elif RCirc == True:
        length4index = len(ABCDmatrices)
        cascTF = matMult_1(ABCDmatrices)

### Compute outputs for circuits including capacitors and/or inductors
if LR == True or RC == True or LCR == True or LC == True:
    for a in range(len(cascTFmatrices)):
        ZinV.append((((termsVector[0][4])*(cascTFmatrices[a][0][0][0]))+cascTFmatrices[a][0][0][1])/(((termsVector[0][4])*(cascTFmatrices[a][0][1][0]))+cascTFmatrices[a][0][1][1]))
        VinV.append(termsVector[0][0]*(ZinV[a]/(ZinV[a]+termsVector[0][2])))
        IinV.append(VinV[a]/ZinV[a])
        PinV.append(VinV[a]*np.conj(IinV[a]))
        inputMatV.append([[VinV[a]],[IinV[a]]])
        AvHelp = 1/((cascTFmatrices[a][0][0][0])+(cascTFmatrices[a][0][0][1])*(1/termsVector[0][4]))
        AiHelp = 1/((cascTFmatrices[a][0][1][1])+(cascTFmatrices[a][0][1][0])*(termsVector[0][4]))
        VoutV.append(VinV[a]*AvHelp)
        IoutV.append(IinV[a]*AiHelp)
        PoutV.append(VoutV[a]*np.conj(IoutV[a]))
        ZoutV.append((((termsVector[0][2])*(cascTFmatrices[a][0][1][1]))+cascTFmatrices[a][0][0][1])/(((termsVector[0][2])*(cascTFmatrices[a][0][1][0]))+cascTFmatrices[a][0][0][0]))
        AvV.append(VoutV[a]/VinV[a])
        AiV.append(IoutV[a]/IinV[a])
        outV.append([VinV[a],IinV[a],VoutV[a],IoutV[a],ZinV[a],ZoutV[a],PinV[a],PoutV[a],AvV[a],AiV[a]])

### Compute outputs for purely resistive circuits
elif RCirc == True:
        invcascTF = np.linalg.inv(cascTF)
        Zin = ((((termsVector[0][4]) * (cascTF[0][0][0])) + cascTF[0][0][1])/(((termsVector[0][4]) * (cascTF[0][1][0])) + cascTF[0][1][1]))
        Vin = (termsVector[0][0] * (Zin / (Zin + termsVector[0][2])))
        Iin = (Vin / Zin)
        Pin = (Vin * Iin)
        inputMat = [[Vin], [Iin]]
        outputM = (np.dot(invcascTF, inputMat))
        Vout = (outputM[0][0])
        Iout = (outputM[0][1])
        Pout = (Vout * Iout)
        Zout = (((termsVector[0][2]) * (cascTF[0][1][1])) + cascTF[0][0][1]) / (((termsVector[0][2]) * (cascTF[0][1][0])) + cascTF[0][0][0])
        Av = (Vout / Vin)
        Ai = ((Iout / Iin))
        outqs.append([complex(Vin, 0), complex(*(Vout, 0)), complex(*(Iin, 0)), complex(*(Iout, 0)), complex(Pin, 0),complex(Zout, 0), complex(*(Pout, 0)), complex(Zin, 0), complex(*(Av, 0)), complex(*(Ai, 0))])

for i in outputVector:
    outputVector = [i.strip(' ') for i in outputVector] # Strip output vector strings of any leading or trailing whitespace
                                                        # to allow for easy reading in write out section

################################################################################################################################################################################################################################################################
######################################## WRITE TO OUTPUT FILE ##################################################################################################################################################################################################
################################################################################################################################################################################################################################################################

### Time response write-out block
if sys.argv[1] == "-i": # Check if time response was called in command line
    ### For purely resistive circuits
    if RCirc == True:
        AvTR = []
        for i in range(len(freqRangeIFFT)):
            AvTR.append(outqs[0][8]) # Av complex number, same for all frequencies due to purely resistive circuit
    ### Otherwise
    elif RC == True or LC == True or LCR == True or LR == True:
        AvV = []
        m=0
        for a in range(len(cascTFmatrices)): # Calculate voltage gain
            AvV.append(1/((cascTFmatrices[a][0][0][0]) + (cascTFmatrices[a][0][0][1]) * (1 / termsVector[0][4])))
        for i in range(len(freqRangeIFFT)):
            if freqRangeIFFT[i] < 0:
                m=m+1
                AvTemp = np.conj(AvV[i]) # For negative frequencies, take complex conjugate of Av
                AvV[i] = AvTemp # Set AvV element with correct value
            else:
                pass
        AvNeg=AvV[0:m] # Extract all the negative frequency terms
        del AvV[0:m] # Delete these terms from AvV
        AvTR = AvV+AvNeg # Concatenate to get correct order (Least positive - Most positive - Most negative - Least negative)

    TimeResp = np.fft.ifft(AvTR) # Calculate time response
    W = open(sys.argv[2] + '.out', 'w+') # Open output file
    ### Write out
    W.write('TIME'.rjust(6))
    W.write(('Vout(t)/Vin(t)').rjust(22))
    W.write('\n')
    W.write('Secs'.rjust(6))
    W.write('L'.rjust(17))
    timestep = 1/(max(freqRangeIFFT)*2) # Calculation of time step
    for i in range(len(freqRangeIFFT)):
        W.write('\n')
        timePoint = timestep*i
        W.write((str.format('{0:.3e}', timePoint)).rjust(9))         # Writing out time steps with
        W.write((str.format('{0:.3e}', TimeResp[i].real)).rjust(11)) # time response
        W.write('+j')
        W.write((str.format('{0:.3e}', TimeResp[i].imag)).rjust(10)) # Get correct complex format
    W.write('\n')
    W.close()

### Frequency response write-out block
else:
    try: # If called from command line
     W= open(sys.argv[2],"w+") # Open write out file
    except: # If called from development environment
     W=open('My_Circuit.dat',"w+") # Open write out file

    if LCR == True or LC == True or RC == True or LR == True:
      frequencyOFInt = []
      if LFREQ == True:
          frequencyOFInt = Lfrequency # Logarithmic frequency sweep
      elif FREQ == True:
          frequencyOFInt = frequency # Normal frequency sweep

    ### Write out
    W.write('Freq'.rjust(10))
    for output in outputVector:
        # Loop through outputVector to find what needs outputting and in what order
        # Use the if loop to pick correct output value
        if output == outputVector[0]: # If first output quantity, need smaller gap
            space = 14
        else:
            space = 23

        ### Write first line containing output title
        if output == 'Vin V':
            W.write('Vin'.rjust(space))

        elif output == 'Vin dBV':
            W.write('Vin'.rjust(space))

        elif output == 'Vout V':
            W.write('Vout'.rjust(space))

        elif output == 'Vout dBV':
            W.write('Vout'.rjust(space))

        elif output == 'Iin A':
            W.write('Iin'.rjust(space))

        elif output == 'Iin dBA':
            W.write('Iin'.rjust(space))

        elif output == 'Iout A':
            W.write('Iout'.rjust(space))

        elif output == 'Iout dBA':
            W.write('Iout'.rjust(space))

        elif output == 'Pin W':
            W.write('Pin'.rjust(space))

        elif output == 'Pin dBW':
            W.write('Pin'.rjust(space))

        elif output == 'Pout W':
            W.write('Pout'.rjust(space))

        elif output == 'Pout dBW':
            W.write('Pout'.rjust(space))

        elif output == 'Zin Ohms':
            W.write('Zin'.rjust(space))

        elif output == 'Zout Ohms':
            W.write('Zout'.rjust(space))

        elif output == 'Av':
            W.write('Av'.rjust(space))

        elif output == 'Av dB':
            W.write('Av'.rjust(space))

        elif output == 'Ai':
            W.write('Ai'.rjust(space))

        elif output == 'Ai dB':
            W.write('Ai'.rjust(space))
        else:
            pass

        if output == outputVector[len(outputVector)-1]:
            W.write('          ') # Add in whitespace found in the .out files
        else:
            pass
    W.write('\n')

    W.write('Hz'.rjust(10))
    for output in outputVector:
    # Loop through outputVector to find desired output quantities and their order
          if output == outputVector[0]: # If first output quantity, need smaller gap
              space = 14
          else:
              space = 23

          ### Write first line containing output units
          if output == 'Vin V':
              W.write('V'.rjust(space))

          elif output == 'Vin dBV':
              W.write('dBV'.rjust(space))

          elif output == 'Vout V':
              W.write('V'.rjust(space))

          elif output == 'Vout dBV':
              W.write('dBV'.rjust(space))

          elif output == 'Iin A':
              W.write('A'.rjust(space))

          elif output == 'Iin dBA':
              W.write('dBA'.rjust(space))

          elif output == 'Iout A':
              W.write('A'.rjust(space))

          elif output == 'Iout dBA':
              W.write('dBA'.rjust(space))

          elif output == 'Pin W':
              W.write('W'.rjust(space))

          elif output == 'Pin dBW':
              W.write('dBW'.rjust(space))

          elif output == 'Pout W':
              W.write('W'.rjust(space))

          elif output == 'Pout dBW':
              W.write('dBW'.rjust(space))

          elif output == 'Zin Ohms':
              W.write('Ohms'.rjust(space))

          elif output == 'Zout Ohms':
              W.write('Ohms'.rjust(space))

          elif output == 'Av':
              W.write('L'.rjust(space))

          elif output == 'Av dB':
              W.write('dB'.rjust(space))

          elif output == 'Ai':
              W.write('L'.rjust(space))

          elif output == 'Ai dB':
              W.write('dB'.rjust(space))
          else:
              pass

          if output == outputVector[len(outputVector) - 1]:
              W.write('          ') # Add additional whitespace needed
          else:
              pass
    W.write('\n')

    if LCR == True or LC == True or RC == True or LR == True: # Only use this write out block if not solely resistor based circuit
        for i in range(len(frequencyOFInt)):
                indextoTerminate = len(outputVector)-1
                ### Write out data
                W.write((str.format('{0:.3e}', frequencyOFInt[i])).rjust(10))
                for item in outputVector:
                # Write out output quantity in correct format
                # calculating dB values if necessary
                   if item  == 'Vin V':
                       W.write((str.format('{0:.3e}', VinV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', VinV[i].imag)).rjust(10))
                   elif item == 'Vin dBV':
                       VinVdB = cmath.polar(VinV[i])
                       VinMod = 20*math.log10(VinVdB[0])
                       VinArg = VinVdB[1]
                       W.write((str.format('{0:.3e}', VinMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', VinArg)).rjust(10))
                   elif item == 'Vout V':
                       W.write((str.format('{0:.3e}', VoutV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', VoutV[i].imag)).rjust(10))
                   elif item == 'Vout dBV':
                       VoutdB = cmath.polar(VoutV[i])
                       VoutMod = 20 * math.log10(VoutdB[0])
                       VoutArg = VoutdB[1]
                       W.write((str.format('{0:.3e}', VoutMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', VoutArg)).rjust(10))
                   elif item == 'Iin A':
                       W.write((str.format('{0:.3e}', IinV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', IinV[i].imag)).rjust(10))
                   elif item == 'Iin dBA':
                       IindB = cmath.polar(IinV[i])
                       IinMod = 20 * math.log10(IindB[0])
                       IinArg = IindB[1]
                       W.write((str.format('{0:.3e}', IinMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', IinArg)).rjust(10))
                   elif item == 'Iout A':
                       W.write((str.format('{0:.3e}', IoutV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', IoutV[i].imag)).rjust(10))
                   elif item == 'Iout dBA':
                       IoutdB = cmath.polar(IoutV[i])
                       IoutMod = 20 * math.log10(IoutdB[0])
                       IoutArg = IoutdB[1]
                       W.write((str.format('{0:.3e}', IoutMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', IoutArg)).rjust(10))
                   elif item == 'Pin W':
                       W.write((str.format('{0:.3e}', PinV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', PinV[i].imag)).rjust(10))
                   elif item == 'Pin dBW':
                       PindB = cmath.polar(PinV[i])
                       PinMod = 10 * math.log10(PindB[0])
                       PinArg = PindB[1]
                       W.write((str.format('{0:.3e}', PinMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', PinArg)).rjust(10))
                   elif item == 'Pout W':
                       W.write((str.format('{0:.3e}', PoutV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', PoutV[i].imag)).rjust(10))
                   elif item == 'Pout dBW':
                       PoutdB = cmath.polar(PoutV[i])
                       PoutMod = 10*math.log10(PoutdBdB[0])
                       PoutArg = PoutdB[1]
                       W.write((str.format('{0:.3e}', PoutMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', PoutArg)).rjust(10))
                   elif item == 'Zin Ohms':
                       W.write((str.format('{0:.3e}', ZinV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', ZinV[i].imag)).rjust(10))
                   elif item == 'Zout Ohms':
                       W.write((str.format('{0:.3e}', ZoutV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', ZoutV[i].imag)).rjust(10))
                   elif item == 'Av':
                       W.write((str.format('{0:.3e}', AvV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', AvV[i].imag)).rjust(10))
                   elif item == 'Av dB':
                       AvdB = cmath.polar(AvV[i])
                       AvMod = 20 * math.log10(AvdB[0])
                       AvArg = AvdB[1]
                       W.write((str.format('{0:.3e}', AvMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', AvArg)).rjust(10))
                   elif item == 'Ai':
                       W.write((str.format('{0:.3e}', AiV[i].real)).rjust(11))
                       W.write('+j')
                       W.write((str.format('{0:.3e}', AiV[i].imag)).rjust(10))
                   elif item == 'Ai dB':
                       AidB = cmath.polar(AiV[i])
                       AiMod = 20*math.log10(AidB[0])
                       AiArg = AidB[1]
                       W.write((str.format('{0:.3e}', AiMod)).rjust(11))
                       W.write('/_')
                       W.write((str.format('{0:.3e}', AiArg)).rjust(10))
                   else:
                       pass

                   if item == outputVector[indextoTerminate]:
                       W.write('\n') # New line if last output quantity on line

    elif RCirc == True: # Only use this is purely resistive circuit
        indextoTerminate = len(outputVector) - 1
        for i in range(len(frequency)):
            W.write((str.format('{0:.3e}', frequency[i])).rjust(10))
            for item in outputVector:
                if item == 'Vin V':
                    W.write((str.format('{0:.3e}', outqs[0][0].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][0].imag)).rjust(10))
                elif item == 'Vin dBV':
                    VinRdB = cmath.polar(outqs[0][0])
                    VinRdBMod = 20 * math.log10(VinRdB[0])
                    VinRdBArg = VinRdB[1]
                    W.write((str.format('{0:.3e}', VinRdBMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', VinRdBArg)).rjust(10))
                elif item == 'Vout V':
                    W.write((str.format('{0:.3e}', outqs[0][1].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][1].imag)).rjust(10))
                elif item == 'Vout dBV':
                    VoutRdB = cmath.polar(outqs[0][11])
                    VoutRMod = 20 * math.log10(VoutRdB[0])
                    VoutRArg = VoutRdB[1]
                    W.write((str.format('{0:.3e}', VoutRMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', VoutRArg)).rjust(10))
                elif item == 'Iin A':
                    W.write((str.format('{0:.3e}', outqs[0][2].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][2].imag)).rjust(10))
                elif item == 'Iin dBA':
                    IindB = cmath.polar(IinV[i])
                    IinMod = 20 * math.log10(IindB[0])
                    IinArg = IindB[1]
                    W.write((str.format('{0:.3e}', IinMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', IinArg)).rjust(10))
                elif item == 'Iout A':
                    W.write((str.format('{0:.3e}', outqs[0][3].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][3].imag)).rjust(10))
                elif item == 'Iout dBA':
                    IoutdB = cmath.polar(IoutV[i])
                    IoutMod = 20 * math.log10(IoutdB[0])
                    IoutArg = IoutdB[1]
                    W.write((str.format('{0:.3e}', IoutMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', IoutArg)).rjust(10))
                elif item == 'Pin W':
                    W.write((str.format('{0:.3e}', outqs[0][4].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][4].imag)).rjust(10))
                elif item == 'Pin dBW':
                    PindB = cmath.polar(PinV[i])
                    PinMod = 10 * math.log10(PindB[0])
                    PinArg = PindB[1]
                    W.write((str.format('{0:.3e}', PinMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', PinArg)).rjust(10))
                elif item == 'Pout W':
                    W.write((str.format('{0:.3e}', outqs[0][6].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][6].imag)).rjust(10))
                elif item == 'Pout dBW':
                    PoutdB = cmath.polar(PoutV[i])
                    PoutMod = 10*math.log10(PoutdBdB[0])
                    PoutArg = PoutdB[1]
                    W.write((str.format('{0:.3e}', PoutMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', PoutArg)).rjust(10))
                elif item == 'Zin Ohms':
                    W.write((str.format('{0:.3e}', outqs[0][7].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][7].imag)).rjust(10))
                elif item == 'Zout Ohms':
                    W.write((str.format('{0:.3e}', outqs[0][5].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][5].imag)).rjust(10))
                elif item == 'Av':
                    W.write((str.format('{0:.3e}', outqs[0][8].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][8].imag)).rjust(10))
                elif item == 'Av dB':
                    AvdB = cmath.polar(AvV[i])
                    AvMod = 20 * math.log10(AvdB[0])
                    AvArg = AvdB[1]
                    W.write((str.format('{0:.3e}', AvMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', AvArg)).rjust(10))
                elif item == 'Ai':
                    W.write((str.format('{0:.3e}', outqs[0][9].real)).rjust(11))
                    W.write('+j')
                    W.write((str.format('{0:.3e}', outqs[0][9].imag)).rjust(10))
                elif item == 'Ai dB':
                    AiRdB = cmath.polar(outqs[0])
                    AiMod = 20*math.log10(AidB[0])
                    AiArg = AidB[1]
                    W.write((str.format('{0:.3e}', AiMod)).rjust(11))
                    W.write('/_')
                    W.write((str.format('{0:.3e}', AiArg)).rjust(10))
                else:
                    pass

                if item == outputVector[indextoTerminate]:
                       W.write('\n') # New line if last output quantity on line
    W.close()

