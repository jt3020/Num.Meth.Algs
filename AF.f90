SUBROUTINE AFT2D (NN,X,Y,NB,MA,MB,NE,ME)
! ADF MESHING ON PLANAR DOMAINS WITH ELEMENT SIZE BASED ON BOUNDARY SEGMENTS
! INPUT: COORDINATES OF NODAL POINTS, {X(I), Y(I), I=1,NN}
! BOUNDARY SEGMENTS, {MA(I),MB(I),I=1,NB}
! OUTPUT: TRIANGULAR ELEMENTS, {ME(I), I=1,3*NE}
! WORKING ARRAY: {XP(I),YP(I),DP(I), I=1,5000} MID-POINTS AND LENGTHS OF SEGMENTS
! {MA(I),MB(I),I=1,1000} CANDIDATE SEGMENTS IN ELEMENT CONSTRUCTION
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 DIMENSION X(*),Y(*),ME(*),MA(*),MB(*)
DIMENSION MS(10),MT(10),XP(10),YP(10),DP(10)
! COMPUTE THE LENGTHS OF THE BOUNDARY SEGMENTS AND THEIR MID-POINTS 
 NP=NB
 DO I=1,NB
IA=MA(I)
IB=MB(I)
 XP(I)=(X(IA)+X(IB))/2
 YP(I)=(Y(IA)+Y(IB))/2
 DP(I)=(X(IB)-X(IA))**2+(Y(IB)-Y(IA))**2
ENDDO

! Write(*,*) "---"
! Write(*,*) MA(1), MA(2), MA(3), MA(4)
! Write(*,*) MB(1), MB(2), MB(3), MB(4)
! Write(*,*) ME
! Write(*,*) MS
! Write(*,*) MT
! Write(*,*) IA
! Write(*,*) IB
! Write(*,*) XP 
! Write(*,*) YP 
! Write(*,*) DP



! PREPARATION WORKS FOR THE BASE SEGMENT, J1-J2 = LAST SEGMENT ON THE FRONT
 5 J3=0
 J1=MA(NB)
 J2=MB(NB)
 Write(*,*) "5", NB, "->", NB-1
 NB=NB-1
 X1=X(J1)
 Y1=Y(J1)
 X2=X(J2)
 Y2=Y(J2)
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

! Write(*,*) "---"
! Write(*,*) A
! Write(*,*) B
! Write(*,*) C
! Write(*,*) DD
! Write(*,*) IA, IB 
! Write(*,*) J1, J2, J3 
! Write(*,*) TOR, RR 
! Write(*,*) XC, YC 
! Write(*,*) X1, Y1 
! Write(*,*) X2, Y2 
! Write(*,*) XM, YM

! Checked against python up to here


! FILTER OFF SEGMENTS TOO FAR AWAY FROM THE BASE SEGMENT
 9 NS=0
 DO 11 I=1,NB
 IA=MA(I)
 IB=MB(I)
!  Write(*,*) 'DPL #',I,NB,':',DPL(X(IA),Y(IA),X(IB),Y(IB),XC,YC)
!  Write(*,*) X(IA),Y(IA)
!  Write(*,*) X(IB),Y(IB)
!  Write(*,*) XC,YC
 IF (DPL(X(IA),Y(IA),X(IB),Y(IB),XC,YC).GT.RR) GOTO 11
 NS=NS+1
 MS(NS)=IA
 MT(NS)=IB
 11 CONTINUE

! Write(*,*) "---" 
! Write(*,*) IA, IB 
! Write(*,*) MA(:8)
! Write(*,*) MB(:8)
! Write(*,*) 
! Write(*,*) MS 
! Write(*,*) MT 
! Write(*,*) NS
! Write(*,*) DPL(X(IA),Y(IA),X(IB),Y(IB),XC,YC)
! Checked against python up to here - agrees only for first loop

 


! DETERMINE CANDIDATE NODES ON THE GENERATION FRONT
 DO 22 I=1,NS
 J=MS(I)
 P=X(J)
 Q=Y(J)
 IF ((P-XC)**2+(Q-YC)**2.GT.RR.OR.A*P+B*Q.LT.C) GOTO 22
 CALL CHKINT (J1,J2,J,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y,*22)
!  Write(*,*) "In", RR
 CALL CIRCLE (X1,Y1,X2,Y2,P,Q,XC,YC,RR)
!  Write(*,*) "Out", RR
 J3=J
 22 CONTINUE

! Write(*,*) "---"
! Write(*,*) MS 
! Write(*,*) X(:4)
! Write(*,*) Y(:4)
! Write(*,*) 
! Write(*,*) J, P, Q 
! Write(*,*) J3
! Write(*,*) NS

! First loop checked up to here

 IF (J3.EQ.0) THEN
 H=SQRT(RR-TOR-DD/4)
 R=SQRT(RR-TOR)
 AREA=SQRT(DD)*(R+H)
 ALPHA=AREA/((R+H)**2+0.75*DD)
 ELSE
 AREA=A*X(J3)+B*Y(J3)+X1*Y2-X2*Y1
 S=DD+(X(J3)-X1)**2+(Y(J3)-Y1)**2+(X(J3)-X2)**2+(Y(J3)-Y2)**2
 ALPHA=SQRT(12.0)*AREA/S
 ENDIF

! Write(*,*) "---"
! Write(*,*) Area
! Write(*,*) S
! Write(*,*) Alpha
! Write(*,*) 
! Write(*,*) J3 
! Write(*,*) A, X(1)
! Write(*,*) B, Y(1)
! Write(*,*) X1, Y1 
! Write(*,*) X2, Y2

! CREATE INTERIOR NODES, CHECK THEIR QUALITIES AND COMPARE WITH FRONTAL NODE J3
 XX=XM+A/2
 YY=YM+B/2
 S1=0
 S2=0
 DO 44 I=1,NP
 S=(XP(I)-XX)**2+(YP(I)-YY)**2+TOR
 S1=S1+DP(I)/S
 S2=S2+1/S !44
 44 CONTINUE
 F=SQRT(0.75*S1/(S2*DD))
 F1=F
 DO 111 I=1,5
 F1=(2*F1**3+3*F)/(3*F1*F1+2.25) !111
 111 CONTINUE
 S=F*DD/AREA
 IF (S.GT.1) S=1/S
 BETA=S*(2-S)*ALPHA
 T=1/ALPHA-SQRT(ABS(1/ALPHA**2-1))

! First loop checked against python
! Write(*,*) "---"
! Write(*,*) XX, YY 
! Write(*,*) S1, S2 
! Write(*,*) NP 
! Write(*,*) 
! Write(*,*) S 
! Write(*,*) F, F1 
! Write(*,*) Beta, Alpha
! Write(*,*) T 




 DO 66 I=1,9
    S=(11-I)*F1/10
    GAMMA=SQRT(3.0)*S*S*(2-S/F)/(S*S*F+0.75*F)
    IF (GAMMA.LT.BETA) Then 

       
        Write(*,*) "GOTO 1"
        GOTO 1
    EndIf
    P=XM+A*S
    Q=YM+B*S
    IF ((P-XC)**2+(Q-YC)**2.GT.RR) Then 
        Write(*,*) "GOTO 66"
        GOTO 66
    EndIf
    CALL CHKINT (J1,J2,0,X1,Y1,X2,Y2,P,Q,NS,MS,MT,X,Y,*66)
    D=(X(MT(1))-X(MS(1)))**2+(Y(MT(1))-Y(MS(1)))**2
    H=DPL(X(MS(1)),Y(MS(1)),X(MT(1)),Y(MT(1)),P,Q)
        DO 99 J=2,NS
            S=DPL(X(MS(J)),Y(MS(J)),X(MT(J)),Y(MT(J)),P,Q)
            IF (S.GE.H) Then 
                Write(*,*) "GOTO 99"
                GOTO 99
            EndIf
            H=S
            D=(X(MT(J))-X(MS(J)))**2+(Y(MT(J))-Y(MS(J)))**2
            99 CONTINUE
    IF (H.GT.D*T**2) Then 
        Write(*,*) "GOTO 3"
        GOTO 3
    EndIf
    66 CONTINUE

 1 II=3*NE

! Checked against python up to here
! Write(*,*) "---"
! Write(*,*) S, F1
! Write(*,*) Beta, Gamma 
! Write(*,*) P, Q 
! Write(*,*) D, H 
! Write(*,*) A, B


! IF NO NODE CAN BE FOUND TO FORM A VALID ELEMENT WITH THE BASE SEGMENT, 
! ENLARGE THE SEARCH RADIUS
 IF (J3.NE.0) GOTO 2
 IF (RR.GT.100*DD) THEN
 WRITE (*,*) '*** Mesh generation failed! ***'
 RETURN
 ENDIF
 XC=XC+XC-XM
 YC=YC+YC-YM
 RR=(XC-X1)**2+(YC-Y1)**2+TOR
 GOTO 9
! NODE J3 IS FOUND TO FORM VALID ELEMENT WITH BASE SEGMENT J1-J2
! UPDATE GENERATION FRONT WITH FRONTAL NODE J3
 2 NE=NE+1
 ME(II+1)=J1
 ME(II+2)=J2
 ME(II+3)=J3

! !Debug here
! Write(*,*) "---"
! Write(*,*) J1, J2, J3 
! Write(*,*) II 
! Write(*,*) NE 
! Write(*,*) ME(:6)


 DO 77 I=1,NB
 IF (MA(I).NE.J3.OR.MB(I).NE.J1) GOTO 77
 MA(I)=MA(NB)
 MB(I)=MB(NB)
 Write(*,*) "77", NB, "->", NB-1
 NB=NB-1
 GOTO 7
 77 CONTINUE


Write(*,*) "7", NB, "->", NB+1
 NB=NB+1
 MA(NB)=J1
 MB(NB)=J3
 7 CONTINUE

!Debug matches here
! Write(*,*) "---"
! Write(*,*) MA(:8)
! Write(*,*) MB(:8) 
! Write(*,*) NB 

 DO 88 I=1,NB
    IF (MA(I).NE.J2.OR.MB(I).NE.J3) GOTO 88
    IF (NB.EQ.1) Then 
        Write(*,*) "MESH GENERATED"
        RETURN
    EndIf
    MA(I)=MA(NB)
    MB(I)=MB(NB)
    Write(*,*) "88", NB, "->", NB-1
    NB=NB-1
    GOTO 5
    88 CONTINUE
Write(*,*) "89", NB, "->", NB+1
 NB=NB+1
 
 MA(NB)=J3
 MB(NB)=J2
 GOTO 5

 Write(*,*) "NEW NODE", NN+1

! INTERIOR NODE NN CREATED, UPDATE GENERATION FRONT WITH INTERIOR NODE NN
 3 NN=NN+1
 X(NN)=P
 Y(NN)=Q
 II=3*NE
 NE=NE+1
 ME(II+1)=J1
 ME(II+2)=J2
 ME(II+3)=NN
 NB=NB+1
 MA(NB)=J1
 MB(NB)=NN
 NB=NB+1
 MA(NB)=NN
 MB(NB)=J2
 GOTO 5
 END


! CALCULATE THE DISTANCE BETWEEN POINT (X3,Y3) TO LINE SEGMENT (X1,Y1)-(X2,Y2)
 FUNCTION DPL(X1,Y1,X2,Y2,X3,Y3)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 R=(X2-X1)**2+(Y2-Y1)**2
 S=(X2-X1)*(X3-X1)+(Y2-Y1)*(Y3-Y1)
 T=(X3-X1)**2+(Y3-Y1)**2
 DPL=T-S*S/R
 IF (S.GT.R) DPL=(X3-X2)**2+(Y3-Y2)**2
 IF (S.LT.0) DPL=T
 END


! CALCULATE THE CIRCUMCIRCLE OF TRIANGLE (X1,Y1),(X2,Y2),(P,Q)
 SUBROUTINE CIRCLE (X1,Y1,X2,Y2,P,Q,XC,YC,RR)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
 END


 SUBROUTINE CHKINT (J1,J2,J,X1,Y1,X2,Y2,P,Q,NB,MA,MB,X,Y,*)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 DIMENSION MA(*),MB(*),X(*),Y(*)
DATA TOL/0.000001/
! Check if there are any intersections between line segment (P,Q)-(X1,Y1)
! and the non-Delaunay segments MA(i)-MB(i), i=1,NB
 C1=Q-Y1
 C2=P-X1
 C=Q*X1-P*Y1
 CC=C1*C1+C2*C2
 TOR=-TOL*CC*CC
 DO 11 I=1,NB
IA=MA(I)
IB=MB(I)
IF (J.EQ.IA.OR.J.EQ.IB.OR.J1.EQ.IA.OR.J1.EQ.IB) GOTO 11
 XA=X(IA)
 YA=Y(IA)
 XB=X(IB)
 YB=Y(IB)
 IF ((C2*YA-C1*XA+C)*(C2*YB-C1*XB+C).GT.TOR) GOTO 11
 H1=YB-YA
 H2=XB-XA
 H=XA*YB-XB*YA
 IF ((H2*Y1-H1*X1+H)*(H2*Q-H1*P+H).LT.TOR) RETURN 1
 11 CONTINUE
! Check if there are any intersections between line segment (P,Q)-(X2,Y2)
! and the non-Delaunay segments MA(i)-MB(i), i=1,NB
 C1=Q-Y2
 C2=P-X2
 C=Q*X2-P*Y2
 CC=C1*C1+C2*C2
 TOR=-TOL*CC*CC
 DO 22 I=1,NB
IA=MA(I)
IB=MB(I)
IF (J.EQ.IA.OR.J.EQ.IB.OR.J2.EQ.IA.OR.J2.EQ.IB) GOTO 22
 XA=X(IA)
 YA=Y(IA)
 XB=X(IB)
 YB=Y(IB)
 IF ((C2*YA-C1*XA+C)*(C2*YB-C1*XB+C).GT.TOR) GOTO 22
 H1=YB-YA
 H2=XB-XA
 H=XA*YB-XB*YA
 IF ((H2*Y2-H1*X2+H)*(H2*Q-H1*P+H).LT.TOR) RETURN 1
 22 CONTINUE
 END

 Program main
    Implicit None
    Real, dimension(1000) :: MS, MT, XP, YP, DP
    Real*8, dimension(1000) :: X, Y
    Integer, dimension(1000) :: ME, MA, MB 
    Integer :: NN, NB, NE, i

    ! NN = 4
    ! X(1) = 0.0; Y(1) = 0.0
    ! X(2) = 1.0; Y(2) = 0.0
    ! X(3) = 1.0; Y(3) = 1.0
    ! X(4) = 0.0; Y(4) = 1.0

    ! NB = 4
    ! MA(1) = 1; MB(1) = 2
    ! MA(2) = 2; MB(2) = 3
    ! MA(3) = 3; MB(3) = 4
    ! MA(4) = 4; MB(4) = 1

    ! Test problem using finer square grid
    NN = 8
    X(1) = 0.0; Y(1) = 0.0
    X(2) = 0.5; Y(2) = 0.0
    X(3) = 1.0; Y(3) = 0.0
    X(4) = 1.0; Y(4) = 0.5
    X(5) = 1.0; Y(5) = 1.0
    X(6) = 0.5; Y(6) = 1.0
    X(7) = 0.0; Y(7) = 1.0
    X(8) = 0.0; Y(8) = 0.5

    NB = 8
    MA(1) = 1; MB(1) = 2
    MA(2) = 2; MB(2) = 3
    MA(3) = 3; MB(3) = 4
    MA(4) = 4; MB(4) = 5
    MA(1) = 5; MB(1) = 6
    MA(2) = 6; MB(2) = 7
    MA(3) = 7; MB(3) = 8
    MA(4) = 8; MB(4) = 1

    NE = 0

    call AFT2D(NN,X,Y,NB,MA,MB,NE,ME)

    ! write(*,*) '--- Results: ---'
    ! write(*,*) 'Number of Nodes:', NN
    ! Do i = 1, NN
    !     write(*,*) 'Node ',i,' Position:',X(i),Y(i)
    ! EndDo
    ! write(*,*) ''
    ! write(*,*) 'Number of Elements:', NE
    ! Do i = 1, 2
    !     write(*,*) 'Triangle ', i,' Nodes:', ME((i-1)*3+1), ME((i-1)*3 + 2), ME((i-1)*3 + 3)
    ! EndDo
 End Program 