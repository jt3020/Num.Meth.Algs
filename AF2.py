import numpy as np
import matplotlib.pyplot as plt

# Advancing Front Python Routine
# Now with python indexing
def AFT2D(NN,X,Y,NB,MA,MB,NE,ME):

    IA = np.zeros(1000,dtype=int)
    IB = np.zeros(1000,dtype=int)
    MS = np.zeros(1000,dtype=int)
    MT = np.zeros(1000,dtype=int)
    XP = np.zeros(1000)
    YP = np.zeros(1000)
    DP = np.zeros(1000)

    # COMPUTE THE LENGTHS OF THE BOUNDARY SEGMENTS AND THEIR MID-POINTS
    NP = NB

    for I in range(1,NB+1):
        IA = MA[I]
        IB = MB[I]
        XP[I] = (X[IA] + X[IB])/2
        YP[I] = (Y[IA] + Y[IB])/2
        DP[I] = (X[IB] - X[IA])**2 + (Y[IB] - Y[IA])**2

    # Logic of code checked up to here
    #Checked against fortran up to here
    # breakpoint()


    # PREPARATION WORKS FOR THE BASE SEGMENT, J1 - J2 = LAST SEGMENT ON THE FRONT
    LoopBreak = False
    while LoopBreak == False:
        # breakpoint()
        ExitLoop = False 
        J3 = 0 # 5
        J1 = MA[NB]
        J2 = MB[NB]
        NB = NB-1
        X1 = X[J1]
        Y1 = Y[J1]
        X2 = X[J2]
        Y2 = Y[J2]
        A = Y1 - Y2
        B = X2 - X1
        DD = A*A + B*B
        TOR = DD/100
        XM = (X1 + X2)/2
        YM = (Y1 + Y2)/2
        RR = 1.25*DD + TOR
        XC = XM + A # ?
        YC = YM + B
        C = X2*Y1 - X1*Y2 + TOR

        # breakpoint()

        # Debugging performed up to here

        # FILTER OFF SEGMENTS TOO FAR AWAY FROM THE BASE SEGMENT
        ValidElement = False
        while ValidElement == False:
            # breakpoint()
            NS = 0 # 9
            for I in range(1,NB+1): ## Loop 11
                IA = MA[I]
                IB = MB[I]
                # print('DPL #',I,NB,':',DPL(X[IA],Y[IA],X[IB],Y[IB],XC,YC))
                if (DPL(X[IA],Y[IA],X[IB],Y[IB],XC,YC) > RR): continue
                NS = NS + 1
                MS[NS] = IA
                MT[NS] = IB
                # 11 Continue
            
            # Debugging performed up to here (DPL not checked)
            # breakpoint()

            # DETERMINE CANDIDATE NODES ON THE GENERATION FRONT
            for I in range(1,NS+1): # 22 Loop
                
                J = MS[I]
                P = X[J]
                Q = Y[J]
                # breakpoint()
                if ((P-XC)**2+(Q-YC)**2 > RR or A*P+B*Q < C): continue
                if CHKINT(J1,J2,J,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y): # Effectively continue if returns False
                    RR, XC, YC = CIRCLE(X1,Y1,X2,Y2,P,Q,XC,YC,RR) # Only executed if returns true
                    J3=J
                # 22 Continue
                # breakpoint()

            # breakpoint()
            
            if (J3 == 0):
                H = np.sqrt(RR-TOR-DD/4)
                R = np.sqrt(RR-TOR)
                AREA = np.sqrt(DD) * (R+H)
                ALPHA = AREA/((R+H)**2 + 0.75*DD)
            else:
                AREA = A*X[J3] + B*Y[J3] + X1*Y2 - X2*Y1
                S = DD + (X[J3] - X1)**2 + (Y[J3] - Y1)**2 + (X[J3] - X2)**2 + (Y[J3] - Y2)**2
                ALPHA = np.sqrt(12.0)*AREA/S
            # breakpoint()


            # CREATE INTERIOR NODES, CHECK THEIR QUALITIES AND COMPARE WITH FRONTAL NODE J3    
            XX = XM + A/2
            YY = YM + B/2
            S1 = 0
            S2 = 0
            for I in range(1,NP+1): # 44 for loop
                S = (XP[I] - XX)**2 + (YP[I] - YY)**2 + TOR
                S1 = S1 + DP[I]/S
                S2 = S2 + 1/S # 44 Continue
            F = np.sqrt(0.75*S1/(S2*DD))
            F1 = F
            for I in range(1,5+1): # 111 for loop
                F1 = (2*F1**3 + 3*F)/(3*F1*F1 + 2.25) # 111
            S = F*DD/AREA
            if (S > 1): S = 1/S
            BETA = S*(2-S)*ALPHA
            T = 1/ALPHA - np.sqrt(abs(1/ALPHA**2 - 1))

            # breakpoint()

            for I in range(1,9+1): # 66 for loop
                S = (11-I) * F1/10
                GAMMA = np.sqrt(3.0)*S*S*(2-S/F)/(S*S*F+0.75*F)
                if (GAMMA < BETA): print('GOTO 1'); break # GOTO 1
                P = XM + A*S
                Q = YM + B*S

                if ((P-XC)**2 + (Q-YC)**2 > RR): print('GOTO 66'); continue # GOTO 66

                if CHKINT(J1,J2,0,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y) == False: print('GOTO 66 chkint'); continue # GOTO 66

                D = (X[MT[1]] - X[MS[1]])**2 + (Y[MT[1]] - Y[MS[1]])**2
                H = DPL(X[MS[1]],Y[MS[1]],X[MT[1]],Y[MT[1]],P,Q)

                for J in range(2,NS+1): # 99 for loop
                    S = DPL(X[MS[J]],Y[MS[J]],X[MT[J]],Y[MT[J]],P,Q)
                    if (S >= H): print('GOTO 99'); continue # GOTO 99
                    H = S
                    D = (X[MT[J]]-X[MS[J]])**2 + (Y[MT[J]]-Y[MS[J]])**2
                # 99 Continue
                if (H > D*T**2):
                    ExitLoop = True 
                    print('GOTO 3'); break  # !!!!!!!!!!!!!!!!!!!1 GOTO 3 - This allows you to skip the GOTO 5 flag 
                # 66 Continue

            # breakpoint()

            if ExitLoop and GAMMA >= BETA: break
            # 1
            II = 3*NE

            # breakpoint()

            # IF NO NODE CAN BE FOUND TO FORM A VALID ELEMENT WITH THE BASE SEGMENT, ENLARGE THE SEARCH RADIUS
            # breakpoint()
            if (J3 == 0): # GOTO 2 if J3 != 0
                
                if (RR > 100*DD):
                    print('*** Mesh generation failed! ***')
                    return
                XC = XC + XC-XM
                YC = YC + YC-YM
                RR = (XC - X1)**2 + (YC - Y1)**2 + TOR
                # GOTO 9 
                # breakpoint()
                ValidElement = False
                
            else:
                # NODE J3 IS FOUND TO FORM VALID ELEMENT WITH BASE SEGMENT J1-J2 - 
                # breakpoint()
                ValidElement = True
            
        # breakpoint()
        if ExitLoop == False:
            # UPDATE GENERATION FRONT WITH FRONTAL NODE J3
            # 2
            NE = NE+1
            ME[II+1] = J1
            ME[II+2] = J2
            ME[II+3] = J3
            FrontUpdate = True
            # breakpoint()
            for I in range(1,NB+1):
                # breakpoint()
                # print('')
                # print(I,'/',NB-1)
                # print(MA[I], '-', J3)
                # print(MB[I], '-', J1)
                if (MA[I] != J3 or MB[I] != J1): continue
                MA[I] = MA[NB]
                MB[I] = MB[NB]
                NB = NB-1
                # GOTO 7
                # breakpoint()
                # print('Break')
                FrontUpdate = False; break
            # 77 Continue
            # breakpoint()
            if FrontUpdate == True:
                NB = NB + 1
                MA[NB] = J1
                MB[NB] = J3
            # breakpoint()
            # Assume no GOTO 5
            LoopBreak = True
            for I in range(1,NB+1): # Position 7, Loop 88
                # breakpoint()
                if (MA[I] != J2 or MB[I] != J3): 
                    continue # GOTO 88
                if (NB == 1): 
                    # breakpoint()
                    print('*** Mesh generation succeeded! ***')
                    return(NN,NE,ME)
                # breakpoint()
                MA[I] = MA[NB]
                MB[I] = MB[NB]
                NB = NB-1
                # breakpoint()
                LoopBreak = False; break  # GOTO 5 ################## 5

            # breakpoint()

            if LoopBreak == False: continue
            # 88 
            NB = NB+1
            MA[NB] = J3
            MB[NB] = J2
            LoopBreak = False
            # print('LOOPBREAK:',LoopBreak)
            # breakpoint()
            continue # GOTO 5 ######################### 5

        # INTERIOR NODE NN CREATED, UPDATE GENERATION FRONT WITH INTERIOR NODE NN
        # 3 
        
        NN=NN+1
        X[NN]=P
        Y[NN]=Q
        II=3*NE
        NE=NE+1
        ME[II+1] = J1
        ME[II+2] = J2
        ME[II+3] = NN
        NB=NB+1
        MA[NB]=J1
        MB[NB]=NN
        NB=NB+1
        MA[NB]=NN
        MB[NB]=J2
        # breakpoint()
        LoopBreak = False; continue # GOTO 5 ################# 5
    # End of Log5 Loop
    
    # End of subroutine
    # breakpoint()
    # print('Finished without Return')

# CALCULATE THE DISTANCE BETWEEN POINT (X3,Y3) TO LINE SEGMENT (X1,Y1)-(X2,Y2)
# Done
# Returns distance DPL
def DPL(X1,Y1,X2,Y2,X3,Y3):
    R = (X2-X1)**2 + (Y2-Y1)**2
    S = (X2-X1)*(X3-X1) + (Y2-Y1)*(Y3-Y1)
    T = (X3-X1)**2 + (Y3-Y1)**2
    
    if (S > R): 
        DPL = (X3-X2)**2 + (Y3-Y2)**2
    elif (S < 0): 
        DPL = T
    else:
        DPL = T-S*S/R
    
    # DPL = np.sqrt(DPL)

    # breakpoint()

    return(DPL)
    # End of Subroutine


# CALCULATE THE CIRCUMCIRCLE OF TRIANGLE (X1,Y1),(X2,Y2),(P,Q)
# Done
# Calculates the coordinates of the centre of the circle and the radius - XC, YC, RR
def CIRCLE(X1,Y1,X2,Y2,P,Q,XC,YC,RR):
    A1=X2-X1
    A2=Y2-Y1
    B1=P-X1
    B2=Q-Y1
    AA=A1*A1+A2*A2
    BB=B1*B1+B2*B2
    AB=A1*B1+A2*B2
    DET=AA*BB-AB*AB
    C1=0.5*BB*(AA-AB)/DET
    C2=0.5*AA*(BB-AB)/DET
    XX=C1*A1+C2*B1
    YY=C1*A2+C2*B2
    RR=1.000001*(XX*XX+YY*YY)
    XC=X1+XX
    YC=Y1+YY
    return(RR,XC,YC)
    # End of subroutine


# '*' removed from arguments
# Converted to python indexing
# Done - Returns True or False
def CHKINT(J1,J2,J,X1,Y1,X2,Y2,P,Q,NB,MA,MB,X,Y):
    TOL = 0.000001

    Check = True

    # Check if there are any intersections between line segment (P,Q)-(X1,Y1)
    # and the non-Delaunay segments MA(i)-MB(i), i=1,NB
    C1 = Q - Y1
    C2 = P - X1
    C = Q*X1 - P*Y1
    CC = C1*C1 + C2*C2
    TOR = -TOL*CC*CC
    for I in range(1,NB+1): #Loop 11
        IA = MA[I]
        IB = MB[I]
        if (J == IA or J == IB or J1 == IA or J1 == IB): continue # GOTO 11
        XA = X[IA]
        YA = Y[IA]
        XB = X[IB]
        YB = Y[IB]
        if ((C2*YA - C1*XA+C)*(C2*YB - C1*XB+C) > TOR): continue # GOTO 11
        H1 = YB - YA
        H2 = XB - XA
        H = XA*YB - XB*YA
        if ((H2*Y1 - H1*X1+H)*(H2*Q - H1*P+H) < TOR): Check = False # RETURN 1 - Returns a fail condition
    # 11 Continue

    # Check if there are any intersections between line segment (P,Q)-(X2,Y2)
    # and the non-Delaunay segments MA(i)-MB(i), i=1,NB    
    C1 = Q - Y2
    C2 = P - X2
    C = Q*X2 - P*Y2
    CC = C1*C1 + C2*C2
    TOR = -TOL*CC*CC
    for I in range(NB): # Loop 22
        IA = MA[I]
        IB = MB[I]
        if (J == IA or J == IB or J2 == IA or J2 == IB): continue # GOTO 22
        XA = X[IA]
        YA = Y[IA]
        XB = X[IB]
        YB = Y[IB]
        if ((C2*YA - C1*XA+C)*(C2*YB - C1*XB+C) > TOR): continue # GOTO 22
        H1 = YB - YA
        H2 = XB - XA
        H = XA*YB - XB*YA
        if ((H2*Y2 - H1*X2+H)*(H2*Q - H1*P+H) < TOR): Check = False # RETURN 1 - Returns a fail condition
    # 22 
    if Check == False:
        return(False)
    else: 
        return(True)
    # End of Subroutine



# NN - Number of nodal points
# X[:], Y[:] - Co-ordinates of nodal points

# NB - Number of boundary segments
# MA[:], MB[:] - Boundary segment nodal IDs

# NE - Number of triangular elements
# ME[:] - Triangular element nodal IDs where ME is of size 3*NE

# XP[:], YP[:], DP[:] - Working arrays containing mid-points and segment lengths - length 5000
# MS[:], MT[:] - Candidate segments used in element construction - length 1000

# Known dimension
MS = np.zeros(1000)
MT = np.zeros(1000)
XP = np.zeros(5000)
YP = np.zeros(5000)
DP = np.zeros(5000)

# Unknown dimension
X = np.zeros(1000)
Y = np.zeros(1000)
ME = np.zeros(1000,dtype=int)
MA = np.zeros(1000,dtype=int)
MB = np.zeros(1000,dtype=int)

# Test problem using square grid
NN = 4
X[1] = 0.0; Y[1] = 0.0
X[2] = 1.0; Y[2] = 0.0
X[3] = 1.0; Y[3] = 1.0
X[4] = 0.0; Y[4] = 1.0

NB = 4
MA[1] = 1; MB[1] = 2
MA[2] = 2; MB[2] = 3
MA[3] = 3; MB[3] = 4
MA[4] = 4; MB[4] = 1

NE = 0

NN, NE, ME = AFT2D(NN,X,Y,NB,MA,MB,NE,ME)

print('--- Results: ---')
print('Number of Nodes:', NN)
for i in range(NN):
    print('Node ',i,' Position:',X[i],Y[i])
print('')
print('Number of Elements:', NE)
for i in range(NE):
    print('Triangle ', i,' Nodes:',ME[i*3 + 1],ME[i*3 + 2],ME[i*3 + 3])






