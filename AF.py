import numpy as np
import matplotlib.pyplot as plt

# Advancing Front Python Routine
# Not yet in python indexing
def AFT2D(NN,X,Y,NB,MA,MB,NE,ME):

    # COMPUTE THE LENGTHS OF THE BOUNDARY SEGMENTS AND THEIR MID-POINTS
    NP = NB
    for I in range(1,NB+1):
        IA = MA[I-1]
        IB = MB[I-1]
        XP[I-1] = (X[IA-1] + X[IB-1])/2
        YP[I-1] = (Y[IA-1] + Y[IB-1])/2
        DP[I-1] = (X[IB-1] - X[IA-1])**2 + (Y[IB-1] - Y[IA-1])**2

    # PREPARATION WORKS FOR THE BASE SEGMENT, J1 - J2 = LAST SEGMENT ON THE FRONT
    LoopBreak = False
    while LoopBreak == False:
        ExitLoop = False 
        J3 = 0 # 5
        J1 = MA[NB-1]
        J2 = MB[NB-1]
        NB = NB-1
        X1 = X[J1-1]
        Y1 = Y[J1-1]
        X2 = X[J2-1]
        Y2 = Y[J2-1]
        A = Y1 - Y2
        B = X2 - X1
        DD = A*A + B*B
        TOR = DD/100
        XM = (X1 + X2)/2
        YM = (Y1 + Y2)/2
        RR = 1.25*DD + TOR
        XC = XM + A
        YC = YM + B
        C = X2*Y1 - X1*Y2 + TOR

        # FILTER OFF SEGMENTS TOO FAR AWAY FROM THE BASE SEGMENT
        ValidElement = False
        while ValidElement == False:
            NS=0 # 9
            for I in range(1,NB+1): ## Loop 11
                IA = MA[I-1]
                IB = MB[I-1]
                if (DPL(X[IA-1],Y[IA-1],X[IB-1],Y[IB-1],XC,YC) > RR): continue
                NS = NS+1
                MS[NS-1] = IA
                MT[NS-1] = IB
                # 11 Continue

            # DETERMINE CANDIDATE NODES ON THE GENERATION FRONT
            for I in range(1,NS+1): # 22 Loop
                J = MS[I-1]
                P = X[J-1]
                Q = Y[J-1]
                if ((P-XC)**2+(Q-YC)**2 > RR or A*P+B*Q < C): continue
                if CHKINT(J1,J2,J,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y): # Effectively continue if returns False
                    CIRCLE(X1,Y1,X2,Y2,P,Q,XC,YC,RR)
                    J3=J
                # 22 Continue
            
            if (J3 == 0):
                H = np.sqrt(RR-TOR-DD/4)
                R = np.sqrt(RR-TOR)
                AREA = np.sqrt(DD) * (R+H)
                ALPHA = AREA/((R+H)**2 + 0.75*DD)
            else:
                AREA = A*X[J3-1] + B*Y[J3-1] + X1*Y2 - X2*Y1
                S = DD + (X[J3-1] - X1)**2 + (Y[J3-1] - Y1)**2 + (X[J3-1] - X2)**2 + (Y[J3-1] - Y2)**2
                ALPHA = np.sqrt(12.0)*AREA/S
            
            # CREATE INTERIOR NODES, CHECK THEIR QUALITIES AND COMPARE WITH FRONTAL NODE J3    
            XX = XM + A/2
            YY = YM + B/2
            S1 = 0
            S2 = 0
            for I in range(1,NP+1): # 44 for loop
                S = (XP[I-1] - XX)**2 + (YP[I-1] - YY)**2 + TOR
                S1 = S1 + DP[I-1]/S
            S2 = S2 + 1/S # 44 Continue
            F = np.sqrt(0.75*S1/(S2*DD))
            F1 = F
            for I in range(1,5+1): # 111 for loop
                F1 = (2*F1**3 + 3*F)/(3*F1*F1 + 2.25) # 111
            S = F*DD/AREA
            if (S > 1): S = 1/S
            BETA = S*(2-S)*ALPHA
            T = 1/ALPHA - np.sqrt(abs(1/ALPHA**2 - 1))

            for I in range(1,9+1): # 66 for loop
                S = (11-I) * F1/10
                GAMMA = np.sqrt(3.0)*S*S*(2-S/F)/(S*S*F+0.75*F)
                if (GAMMA < BETA): break # GOTO 1
                P = XM + A*S
                Q = YM + B*S

                if ((P-XC)**2 + (Q-YC)**2 > RR): continue # GOTO 66

                if CHKINT(J1,J2,0,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y) == False: continue # GOTO 66

                D = (X(MT[0]) - X(MS[0]))**2 + (Y(MT[0]) - Y(MS[0]))**2
                H = DPL(X(MS[0]),Y(MS[0]),X(MT[0]),Y(MT[0]),P,Q)

                for J in range(2,NS+1): # 99 for loop
                    S = DPL(X(MS[J-1]),Y(MS[J-1]),X(MT[J-1]),Y(MT[J-1]),P,Q)
                    if (S >= H): continue # GOTO 99
                    H = S
                    D = (X(MT[J-1])-X(MS[J-1]))**2 + (Y(MT[J-1])-Y(MS[J-1]))**2
                # 99 Continue
                if (H > D*T**2):
                    ExitLoop = True 
                    break  # !!!!!!!!!!!!!!!!!!!1 GOTO 3 - This allows you to skip the GOTO 5 flag 
                # 66 Continue

            if ExitLoop and GAMMA >= BETA: break
            # 1
            II = 3*NE

            # IF NO NODE CAN BE FOUND TO FORM A VALID ELEMENT WITH THE BASE SEGMENT, ENLARGE THE SEARCH RADIUS
            if (J3 == 0): # GOTO 2 if J3 != 0
                
                if (RR > 100*DD):
                    print('*** Mesh generation failed! ***')
                    return
                XC = XC + XC-XM
                YC = YC + YC-YM
                RR = (XC - X1)**2 + (YC - Y1)**2 + TOR
                # GOTO 9 
                ValidElement = False
                
            else:
                # NODE J3 IS FOUND TO FORM VALID ELEMENT WITH BASE SEGMENT J1-J2 - 
                ValidElement = True
            
        if ExitLoop == False:
            # UPDATE GENERATION FRONT WITH FRONTAL NODE J3
            # 2
            NE = NE+1
            ME[II+0] = J1
            ME[II+1] = J2
            ME[II+2] = J3
            FrontUpdate = True
            for I in range(1,NB+1):
                if (MA[I] != J3 or MB[I] != J1): continue
                MA[I-1] = MA[NB-1]
                MB[I-1] = MB[NB-1]
                NB = NB-1
                # GOTO 7
                if I == NB: FrontUpdate = False
            # 77 Continue
            if FrontUpdate == True:
                NB = NB + 1
                MA[NB-1] = J1
                MB[NB-1] = J3
            for I in range(1,NB+1): # Position 7, Loop 88
                if (MA[I-1] != J2 or MB[I-1] != J3): continue # GOTO 88
                if (NB == 1): 
                    print('*** Mesh generation succeeded! ***')
                    return
                MA[I-1] = MA[NB-1]
                MB[I-1] = MB[NB-1]
                NB = NB-1
                LoopBreak = True # GOTO 5 ################## 5

            if LoopBreak: LoopBreak = False; break
            # 88 
            NB = NB+1
            MA[NB-1] = J3
            MB[NB-1] = J2
            break # GOTO 5 ######################### 5

        # INTERIOR NODE NN CREATED, UPDATE GENERATION FRONT WITH INTERIOR NODE NN
        # 3 
        NN=NN+1
        X[NN-1]=P
        Y[NN-1]=Q
        II=3*NE
        NE=NE+1
        ME[II+0] = J1
        ME[II+1] = J2
        ME[II+2] = NN
        NB=NB+1
        MA[NB-1]=J1
        MB[NB-1]=NN
        NB=NB+1
        MA[NB-1]=NN
        MB[NB-1]=J2
        break # GOTO 5 ################# 5
    # End of Log5 Loop
    # End of subroutine

# CALCULATE THE DISTANCE BETWEEN POINT (X3,Y3) TO LINE SEGMENT (X1,Y1)-(X2,Y2)
# Done
# Returns distance DPL
def DPL(X1,Y1,X2,Y2,X3,Y3):
    R = (X2-X1)**2 + (Y2-Y1)**2
    S = (X2-X1)*(X3-X1) + (Y2-Y1)*(Y3-Y1)
    T = (X3-X1)**2 + (Y3-Y1)**2
    DPL = T-S*S/R
    if (S > R): DPL = (X3-X2)**2 + (Y3-Y2)**2
    if (S < 0): DPL = T
    return(DPL)
    # End of Subroutine


# CALCULATE THE CIRCUMCIRCLE OF TRIANGLE (X1,Y1),(X2,Y2),(P,Q)
# Done
# Calculates the coordinates of the centre of the circle and the radius - XC, YC, RR
def CIRCLE(X1,Y1,X2,Y2,P,Q,XC,YC,RR):
    A1 = X2-X1
    A2 = Y2-Y1
    B1 = P-X1
    B2 = Q-Y1
    AA = A1*A1 + A2*A2
    BB = B1*B1 + B2*B2
    AB = A1*B1 + A2*B2
    DET = AA*BB - AB*AB
    C1 = 0.5*BB*(AA-AB)/DET
    C2 = 0.5*AA*(BB-AB)/DET
    XX = C1*A1 + C2*B1
    YY = C1*A2 + C2*B2
    RR = 1.000001*(XX*XX+YY*YY)
    XC = X1 + XX
    YC = Y1 + YY
    # End of subroutine


# '*' removed from arguments
# Converted to python indexing
# Done - Returns True or False
def CHKINT(J1,J2,J,X1,Y1,X2,Y2,P,Q,NB,MA,MB,X,Y):
    TOL = 0.000001

    # Check if there are any intersections between line segment (P,Q)-(X1,Y1)
    # and the non-Delaunay segments MA(i)-MB(i), i=1,NB
    C1 = Q - Y1
    C2 = P - X1
    C = Q*X1 - P*Y1
    CC = C1*C1 + C2*C2
    TOR = -TOL*CC*CC
    for I in range(NB): #Loop 11
        IA = MA[I]
        IB = MB[I]
        if (J == IA or J == IB or J1 == IA or J1 == IB): break # GOTO 11
        XA = X[IA]
        YA = Y[IA]
        XB = X[IB]
        YB = Y[IB]
        if ((C2*YA - C1*XA+C)*(C2*YB - C1*XB+C) > TOR): break # GOTO 11
        H1 = YB - YA
        H2 = XB - XA
        H = XA*YB - XB*YA
        if ((H2*Y1 - H1*X1+H)*(H2*Q - H1*P+H) < TOR): return(False) # RETURN 1 - Returns a fail condition
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
        if (J == IA or J == IB or J2 == IA or J2 == IB): break # GOTO 22
        XA = X[IA]
        YA = Y[IA]
        XB = X[IB]
        YB = Y[IB]
        if ((C2*YA - C1*XA+C)*(C2*YB - C1*XB+C) > TOR): break # GOTO 22
        H1 = YB - YA
        H2 = XB - XA
        H = XA*YB - XB*YA
        if ((H2*Y2 - H1*X2+H)*(H2*Q - H1*P+H) < TOR): return(False) # RETURN 1 - Returns a fail condition
    # 22 
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

MS = np.zeros(1000)
MT = np.zeros(1000)
XP = np.zeros(5000)
YP = np.zeros(5000)
DP = np.zeros(5000)

X = []
Y = []
ME = []
MA = []
MB = []






