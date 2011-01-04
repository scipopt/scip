*NAME:         enigma
*ROWS:         21
*COLUMNS:      100
*INTEGER:      100
*NONZERO:      289
*BEST SOLN:    0.0 (opt)
*LP SOLN:      0.0
*SOURCE:       Harlan Crowder (IBM)         
*              Harlan Crowder (IBM)
*              E. Andrew Boyd (Rice University)
*APPLICATION:  unknown
*COMMENTS:     pure 0/1 IP
*
*      
NAME          ENIGMA
ROWS
 N  OBJECT  
 E  SOS0    
 E  SOS1    
 E  SOS2    
 E  SOS3    
 E  SOS4    
 E  SOS5    
 E  SOS6    
 E  SOS7    
 E  SOS8    
 E  SOS9    
 E  BILANCIO
 E  SOSA    
 E  SOSB    
 E  SOSC    
 E  SOSD    
 E  SOSE    
 E  SOSF    
 E  SOSG    
 E  SOSH    
 E  SOSI    
 E  SOSL    
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    A0        SOS0                 1   SOSA                 1
    A1        OBJECT               1   SOS1                 1
    A1        BILANCIO           202   SOSA                 1
    A2        OBJECT               2   SOS2                 1
    A2        BILANCIO           404   SOSA                 1
    A3        OBJECT               3   SOS3                 1
    A3        BILANCIO           606   SOSA                 1
    A4        OBJECT               4   SOS4                 1
    A4        BILANCIO           808   SOSA                 1
    A5        OBJECT               5   SOS5                 1
    A5        BILANCIO          1010   SOSA                 1
    A6        OBJECT               6   SOS6                 1
    A6        BILANCIO          1212   SOSA                 1
    A7        OBJECT               7   SOS7                 1
    A7        BILANCIO          1414   SOSA                 1
    A8        OBJECT               8   SOS8                 1
    A8        BILANCIO          1616   SOSA                 1
    A9        OBJECT               9   SOS9                 1
    A9        BILANCIO          1818
    B0        SOS0                 1   SOSB                 1
    B1        SOS1                 1   BILANCIO           -79
    B1        SOSB                 1
    B2        SOS2                 1   BILANCIO          -158
    B2        SOSB                 1
    B3        SOS3                 1   BILANCIO          -237
    B3        SOSB                 1
    B4        SOS4                 1   BILANCIO          -316
    B4        SOSB                 1
    B5        SOS5                 1   BILANCIO          -395
    B5        SOSB                 1
    B6        SOS6                 1   BILANCIO          -474
    B6        SOSB                 1
    B7        SOS7                 1   BILANCIO          -553
    B7        SOSB                 1
    B8        SOS8                 1   BILANCIO          -632
    B8        SOSB                 1
    B9        SOS9                 1   BILANCIO          -711
    B9        SOSB                 1
    C0        SOS0                 1   SOSC                 1
    C1        SOS1                 1   BILANCIO        100023
    C1        SOSC                 1
    C2        SOS2                 1   BILANCIO        200046
    C2        SOSC                 1
    C3        SOS3                 1   BILANCIO        300069
    C3        SOSC                 1
    C4        SOS4                 1   BILANCIO        400092
    C4        SOSC                 1
    C5        SOS5                 1   BILANCIO        500115
    C5        SOSC                 1
    C6        SOS6                 1   BILANCIO        600138
    C6        SOSC                 1
    C7        SOS7                 1   BILANCIO        700161
    C7        SOSC                 1
    C8        SOS8                 1   BILANCIO        800184
    C8        SOSC                 1
    C9        SOS9                 1   BILANCIO        900207
    C9        SOSC                 1
    D0        SOS0                 1   SOSD                 1
    D1        SOS1                 1   BILANCIO        -89810
    D1        SOSD                 1
    D2        SOS2                 1   BILANCIO       -179620
    D2        SOSD                 1
    D3        SOS3                 1   BILANCIO       -269430
    D3        SOSD                 1
    D4        SOS4                 1   BILANCIO       -359240
    D4        SOSD                 1
    D5        SOS5                 1   BILANCIO       -449050
    D5        SOSD                 1
    D6        SOS6                 1   BILANCIO       -538860
    D6        SOSD                 1
    D7        SOS7                 1   BILANCIO       -628670
    D7        SOSD                 1
    D8        SOS8                 1   BILANCIO       -718480
    D8        SOSD                 1
    D9        SOS9                 1   BILANCIO       -808290
    D9        SOSD                 1
    E_0       SOS0                 1   SOSE                 1
    E_1       SOS1                 1   BILANCIO         -9980
    E_1       SOSE                 1
    E_2       SOS2                 1   BILANCIO        -19960
    E_2       SOSE                 1
    E_3       SOS3                 1   BILANCIO        -29940
    E_3       SOSE                 1
    E_4       SOS4                 1   BILANCIO        -39920
    E_4       SOSE                 1
    E_5       SOS5                 1   BILANCIO        -49900
    E_5       SOSE                 1
    E_6       SOS6                 1   BILANCIO        -59880
    E_6       SOSE                 1
    E_7       SOS7                 1   BILANCIO        -69860
    E_7       SOSE                 1
    E_8       SOS8                 1   BILANCIO        -79840
    E_8       SOSE                 1
    E_9       SOS9                 1   BILANCIO        -89820
    E_9       SOSE                 1
    F0        SOS0                 1   SOSF                 1
    F1        SOS1                 1   BILANCIO          1000
    F1        SOSF                 1
    F2        SOS2                 1   BILANCIO          2000
    F2        SOSF                 1
    F3        SOS3                 1   BILANCIO          3000
    F3        SOSF                 1
    F4        SOS4                 1   BILANCIO          4000
    F4        SOSF                 1
    F5        SOS5                 1   BILANCIO          5000
    F5        SOSF                 1
    F6        SOS6                 1   BILANCIO          6000
    F6        SOSF                 1
    F7        SOS7                 1   BILANCIO          7000
    F7        SOSF                 1
    F8        SOS8                 1   BILANCIO          8000
    F8        SOSF                 1
    F9        SOS9                 1   BILANCIO          9000
    F9        SOSF                 1
    G0        SOS0                 1   SOSG                 1
    G1        SOS1                 1   BILANCIO           100
    G1        SOSG                 1
    G2        SOS2                 1   BILANCIO           200
    G2        SOSG                 1
    G3        SOS3                 1   BILANCIO           300
    G3        SOSG                 1
    G4        SOS4                 1   BILANCIO           400
    G4        SOSG                 1
    G5        SOS5                 1   BILANCIO           500
    G5        SOSG                 1
    G6        SOS6                 1   BILANCIO           600
    G6        SOSG                 1
    G7        SOS7                 1   BILANCIO           700
    G7        SOSG                 1
    G8        SOS8                 1   BILANCIO           800
    G8        SOSG                 1
    G9        SOS9                 1   BILANCIO           900
    G9        SOSG                 1
    H0        SOS0                 1   SOSH                 1
    H1        SOS1                 1   BILANCIO         10000
    H1        SOSH                 1
    H2        SOS2                 1   BILANCIO         20000
    H2        SOSH                 1
    H3        SOS3                 1   BILANCIO         30000
    H3        SOSH                 1
    H4        SOS4                 1   BILANCIO         40000
    H4        SOSH                 1
    H5        SOS5                 1   BILANCIO         50000
    H5        SOSH                 1
    H6        SOS6                 1   BILANCIO         60000
    H6        SOSH                 1
    H7        SOS7                 1   BILANCIO         70000
    H7        SOSH                 1
    H8        SOS8                 1   BILANCIO         80000
    H8        SOSH                 1
    H9        SOS9                 1   BILANCIO         90000
    H9        SOSH                 1
    I0        SOS0                 1   SOSI                 1
    I1        SOS1                 1   BILANCIO           100
    I1        SOSI                 1
    I2        SOS2                 1   BILANCIO           200
    I2        SOSI                 1
    I3        SOS3                 1   BILANCIO           300
    I3        SOSI                 1
    I4        SOS4                 1   BILANCIO           400
    I4        SOSI                 1
    I5        SOS5                 1   BILANCIO           500
    I5        SOSI                 1
    I6        SOS6                 1   BILANCIO           600
    I6        SOSI                 1
    I7        SOS7                 1   BILANCIO           700
    I7        SOSI                 1
    I8        SOS8                 1   BILANCIO           800
    I8        SOSI                 1
    I9        SOS9                 1   BILANCIO           900
    I9        SOSI                 1
    L0        SOS0                 1   SOSL                 1
    L1        SOS1                 1   BILANCIO            -1
    L1        SOSL                 1
    L2        SOS2                 1   BILANCIO            -2
    L2        SOSL                 1
    L3        SOS3                 1   BILANCIO            -3
    L3        SOSL                 1
    L4        SOS4                 1   BILANCIO            -4
    L4        SOSL                 1
    L5        SOS5                 1   BILANCIO            -5
    L5        SOSL                 1
    L6        SOS6                 1   BILANCIO            -6
    L6        SOSL                 1
    L7        SOS7                 1   BILANCIO            -7
    L7        SOSL                 1
    L8        SOS8                 1   BILANCIO            -8
    L8        SOSL                 1
    L9        SOS9                 1   BILANCIO            -9
    L9        SOSL                 1
    MARK0001  'MARKER'                 'INTEND'
RHS
    RHS1      SOS0                 1   SOS1                 1
    RHS1      SOS2                 1   SOS3                 1
    RHS1      SOS4                 1   SOS5                 1
    RHS1      SOS6                 1   SOS7                 1
    RHS1      SOS8                 1   SOS9                 1
    RHS1      SOSA                 1   SOSB                 1
    RHS1      SOSC                 1   SOSD                 1
    RHS1      SOSE                 1   SOSF                 1
    RHS1      SOSG                 1   SOSH                 1
    RHS1      SOSI                 1   SOSL                 1
BOUNDS
 UP ONE       A0                   1
 UP ONE       A1                   1
 UP ONE       A2                   1
 UP ONE       A3                   1
 UP ONE       A4                   1
 UP ONE       A5                   1
 UP ONE       A6                   1
 UP ONE       A7                   1
 UP ONE       A8                   1
 UP ONE       A9                   1
 UP ONE       B0                   1
 UP ONE       B1                   1
 UP ONE       B2                   1
 UP ONE       B3                   1
 UP ONE       B4                   1
 UP ONE       B5                   1
 UP ONE       B6                   1
 UP ONE       B7                   1
 UP ONE       B8                   1
 UP ONE       B9                   1
 UP ONE       C0                   1
 UP ONE       C1                   1
 UP ONE       C2                   1
 UP ONE       C3                   1
 UP ONE       C4                   1
 UP ONE       C5                   1
 UP ONE       C6                   1
 UP ONE       C7                   1
 UP ONE       C8                   1
 UP ONE       C9                   1
 UP ONE       D0                   1
 UP ONE       D1                   1
 UP ONE       D2                   1
 UP ONE       D3                   1
 UP ONE       D4                   1
 UP ONE       D5                   1
 UP ONE       D6                   1
 UP ONE       D7                   1
 UP ONE       D8                   1
 UP ONE       D9                   1
 UP ONE       E_0                  1
 UP ONE       E_1                  1
 UP ONE       E_2                  1
 UP ONE       E_3                  1
 UP ONE       E_4                  1
 UP ONE       E_5                  1
 UP ONE       E_6                  1
 UP ONE       E_7                  1
 UP ONE       E_8                  1
 UP ONE       E_9                  1
 UP ONE       F0                   1
 UP ONE       F1                   1
 UP ONE       F2                   1
 UP ONE       F3                   1
 UP ONE       F4                   1
 UP ONE       F5                   1
 UP ONE       F6                   1
 UP ONE       F7                   1
 UP ONE       F8                   1
 UP ONE       F9                   1
 UP ONE       G0                   1
 UP ONE       G1                   1
 UP ONE       G2                   1
 UP ONE       G3                   1
 UP ONE       G4                   1
 UP ONE       G5                   1
 UP ONE       G6                   1
 UP ONE       G7                   1
 UP ONE       G8                   1
 UP ONE       G9                   1
 UP ONE       H0                   1
 UP ONE       H1                   1
 UP ONE       H2                   1
 UP ONE       H3                   1
 UP ONE       H4                   1
 UP ONE       H5                   1
 UP ONE       H6                   1
 UP ONE       H7                   1
 UP ONE       H8                   1
 UP ONE       H9                   1
 UP ONE       I0                   1
 UP ONE       I1                   1
 UP ONE       I2                   1
 UP ONE       I3                   1
 UP ONE       I4                   1
 UP ONE       I5                   1
 UP ONE       I6                   1
 UP ONE       I7                   1
 UP ONE       I8                   1
 UP ONE       I9                   1
 UP ONE       L0                   1
 UP ONE       L1                   1
 UP ONE       L2                   1
 UP ONE       L3                   1
 UP ONE       L4                   1
 UP ONE       L5                   1
 UP ONE       L6                   1
 UP ONE       L7                   1
 UP ONE       L8                   1
 UP ONE       L9                   1
ENDATA
