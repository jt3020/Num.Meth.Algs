import numpy as np
import matplotlib.pyplot as plt

# Advancing Front Python Routine

def AFT2D(NN,X,Y,NB,MA,MB,NE,ME):
    # COMPUTE THE LENGTHS OF THE BOUNDARY SEGMENTS AND THEIR MID-POINTS
    NP = NB
    for I in range(1,NB):
        IA = MA[I]
        IB = MB[I]
        XP[I] = (X[IA]+X[IB])/2
        YP[I] = (Y[IA]+Y[IB])/2
        DP[I] = (X[IB]-X[IA])**2+(Y[IB]-Y[IA])**2

    # PREPARATION WORKS FOR THE BASE SEGMENT, J1-J2 = LAST SEGMENT ON THE FRONT
    J3 = 0
    J1 = MA[NB]
    J2 = MB[NB]
    NB = NB-1
    X1 = X[J1]
    Y1 = Y[J1]
    X2=X[J2]
    Y2=Y[J2]
    A=Y1-Y2
    B=X2-X1
    DD=A*A+B*B
    TOR=DD/100
    XM=(X1+X2)/2
    YM=(Y1+Y2)/2
    RR=1.25*DD+TOR
    XC=XM+A
    YC=YM+B
    C=X2*Y1-X1*Y2+TOR

    # FILTER OFF SEGMENTS TOO FAR AWAY FROM THE BASE SEGMENT
    NS=0
    for I in range(1,NB):
        IA = MA[I]
        IB = MB[I]
        # if (DPL(X[IA],Y[IA],X[IB],Y[IB],XC,YC) > RR): break
        NS = NS+1
        MS[NS] = IA
        MT[NS] = IB

    # DETERMINE CANDIDATE NODES ON THE GENERATION FRONT
    for I in range(1,NS):
        J = MS[I]
        P = X[J]
        Q = Y[J]
        if ((P-XC)**2+(Q-YC)**2 > RR or A*P+B*Q < C): break
        # CHKINT(J1,J2,J,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y,*22)
        # CIRCLE(X1,Y1,X2,Y2,P,Q,XC,YC,RR)
        J3=J
    
    if (J3 == 0):
        H=np.sqrt(RR-TOR-DD/4)
        R=np.sqrt(RR-TOR)
        AREA=np.sqrt(DD)*(R+H)
        ALPHA=AREA/((R+H)**2+0.75*DD)
    else:
        AREA=A*X[J3]+B*Y[J3]+X1*Y2-X2*Y1
        S=DD+(X[J3]-X1)**2+(Y[J3]-Y1)**2+(X[J3]-X2)**2+(Y[J3]-Y2)**2
        ALPHA=np.sqrt(12.0)*AREA/S
    
    # CREATE INTERIOR NODES, CHECK THEIR QUALITIES AND COMPARE WITH FRONTAL NODE J3    
    XX = XM+A/2
    YY = YM+B/2
    S1 = 0
    S2 = 0
    for I in range(1,NP):
        S=(XP[I]-XX)**2+(YP[I]-YY)**2+TOR
        S1=S1+DP[I]/S
    S2=S2+1/S
    F=np.sqrt(0.75*S1/(S2*DD))
    F1=F
    for I in range(1,5): 
        F1=(2*F1**3+3*F)/(3*F1*F1+2.25)
    S=F*DD/AREA
    if (S > 1): S=1/S
    BETA=S*(2-S)*ALPHA
    T=1/ALPHA-np.sqrt(abs(1/ALPHA**2-1))

    for I in range(1,9):
        S=(11-I)*F1/10
        GAMMA=np.sqrt(3.0)*S*S*(2-S/F)/(S*S*F+0.75*F)
        if (GAMMA < BETA): break # GOTO 1
        P=XM+A*S
        Q=YM+B*S
        if ((P-XC)**2+(Q-YC)**2 > RR): break

        # CHKINT(J1,J2,0,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y,*66)

        D = (X(MT[1])-X(MS[1]))**2 + (Y(MT[1])-Y(MS[1]))**2
        # H = DPL(X(MS[1]),Y(MS[1]),X(MT[1]),Y(MT[1]),P,Q)
        for J in range(2,NS):
            # S=DPL(X(MS[J]),Y(MS[J]),X(MT[J]),Y(MT[J]),P,Q)
            if (S >= H): break # GOTO 99
            H=S
            D=(X(MT[J])-X(MS[J]))**2+(Y(MT[J])-Y(MS[J]))**2
        if (H > D*T**2): break  # GOTO 3
    # 1
    II=3*NE

    # IF NO NODE CAN BE FOUND TO FORM A VALID ELEMENT WITH THE BASE SEGMENT, ENLARGE THE SEARCH RADIUS
    if (J3 != 0): pass
        # GOTO 2
        
    if (RR > 100*DD):
        print('*** Mesh generation failed! ***')
        return
    XC = XC+XC-XM
    YC = YC+YC-YM
    RR = (XC-X1)**2 + (YC-Y1)**2+TOR
    # GOTO 9
    # NODE J3 IS FOUND TO FORM VALID ELEMENT WITH BASE SEGMENT J1-J2 UPDATE GENERATION FRONT WITH FRONTAL NODE J3
    # 2
    NE = NE+1
    ME[II+1] = J1
    ME[II+2] = J2
    ME[II+3] = J3
    for I in range(1,NB):
        if (MA(I) != J3 or MB(I) != J1): continue
        MA[I] = MA[NB]
        MB[I] = MB[NB]
        NB = NB-1
        # GOTO 7
    # 77 Continue
    NB = NB 
    MA[NB] = J1
    MB[NB] = J3
    for I in range(1,NB): # Position 7, Loop 88
        if (MA(I) != J2 or MB(I) != J3): break# GOTO 88
        if (NB == 1): return
        MA[I] = MA[NB]
        MB[I] = MB[NB]
        NB = NB-1
        # GOTO 5
    # 88 
    NB = NB+1
    MA[NB] = J3
    MB[NB] = J2
    # GOTO 5
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
    # GOTO 5
    # End of subroutine

# CALCULATE THE DISTANCE BETWEEN POINT (X3,Y3) TO LINE SEGMENT (X1,Y1)-(X2,Y2)
def DPL(X1,Y1,X2,Y2,X3,Y3):
    R=(X2-X1)**2+(Y2-Y1)**2
    S=(X2-X1)*(X3-X1)+(Y2-Y1)*(Y3-Y1)
    T=(X3-X1)**2+(Y3-Y1)**2
    DPL=T-S*S/R
    if (S > R): DPL=(X3-X2)**2+(Y3-Y2)**2
    if (S < 0): DPL=T
    # End of Subroutine


# CALCULATE THE CIRCUMCIRCLE OF TRIANGLE (X1,Y1),(X2,Y2),(P,Q)
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
    # End of subroutine

# Check if there are any intersections between line segment (P,Q)-(X1,Y1)
# and the non-Delaunay segments MA(i)-MB(i), i=1,NB
def CHKINT(J1,J2,J,X1,Y1,X2,Y2,P,Q,NB,MA,MB,X,Y,*):
    TOL=0.000001
    C1=Q-Y1
    C2=P-X1
    C=Q*X1-P*Y1
    CC=C1*C1+C2*C2
    TOR=-TOL*CC*CC
    for I in range(1,NB): #Loop 11
        IA = MA[I]
        IB = MB[I]
        if (J == IA or J == IB or J1 == IA or J1 == IB): break # GOTO 11
        XA = X[IA]
        YA = Y[IA]
        XB = X[IB]
        YB = Y[IB]
        if ((C2*YA-C1*XA+C)*(C2*YB-C1*XB+C) > TOR): break # GOTO 11
        H1 = YB-YA
        H2 = XB-XA
        H = XA*YB-XB*YA
        if ((H2*Y1-H1*X1+H)*(H2*Q-H1*P+H) < TOR): return # RETURN 1
    # 11 Continue
    # Check if there are any intersections between line segment (P,Q)-(X2,Y2)
    # and the non-Delaunay segments MA(i)-MB(i), i=1,NB
        
    C1=Q-Y2
    C2=P-X2
    C=Q*X2-P*Y2
    CC=C1*C1+C2*C2
    TOR=-TOL*CC*CC
    for I in range(1,NB): # Loop 22
        IA = MA[I]
        IB = MB[I]
        if (J == IA or J == IB or J2 == IA or J2 == IB): break # GOTO 22
        XA=X(IA)
        YA=Y(IA)
        XB=X(IB)
        YB=Y(IB)
        if ((C2*YA-C1*XA+C)*(C2*YB-C1*XB+C) > TOR): break # GOTO 22
        H1=YB-YA
        H2=XB-XA
        H=XA*YB-XB*YA
        if ((H2*Y2-H1*X2+H)*(H2*Q-H1*P+H) < TOR): return # RETURN 1
    # 22 
    # End of Subroutine



    




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






