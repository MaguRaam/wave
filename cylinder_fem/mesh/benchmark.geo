//Domain:
m=1.0;
//Geometry:
L=4;
H=2;
 
Point(1) = {0, 0, 0, m};
Point(2) = {L, 0, 0, m};
Point(3) = {L, H, 0, m};
Point(4) = {0, H, 0, m};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//Array of cylinders:

//radius:

//Template cylinder;
R=0.15;
C1=H/4;

Point(5) = {L/4,C1 , 0, m};
Point(6) = {L/4,(C1)+R , 0, m};
Point(7) = {L/4,(C1)-R , 0, m};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 6};

Point(8) = {L/4,2*C1 , 0, m};
Point(9) = {L/4,2*(C1)+R , 0, m};
Point(10) = {L/4,2*(C1)-R , 0, m};

Circle(7) = {10, 8, 9};
Circle(8) = {9, 8, 10};

Point(11) = {L/4,3*C1 , 0, m};
Point(12) = {L/4,3*(C1)+R , 0, m};
Point(13) = {L/4,3*(C1)-R , 0, m};
 
Circle(9) = {13, 11, 12};
Circle(10) = {12, 11, 13};



Point(15) = {L/2,C1 , 0, m};
Point(16) = {L/2,(C1)+R , 0, m};
Point(17) = {L/2,(C1)-R , 0, m};

Circle(15) = {16, 15, 17};
Circle(16) = {17, 15, 16};

Point(18) = {L/2,2*C1 , 0, m};
Point(19) = {L/2,2*(C1)+R , 0, m};
Point(110) = {L/2,2*(C1)-R , 0, m};

Circle(17) = {110, 18, 19};
Circle(18) = {19, 18, 110};

Point(111) = {L/2,3*C1 , 0, m};
Point(112) = {L/2,3*(C1)+R , 0, m};
Point(113) = {L/2,3*(C1)-R , 0, m};
 
Circle(19) = {113, 111, 112};
Circle(110) = {112, 111, 113};


Point(25) = {3*L/4,C1 , 0, m};
Point(26) = {3*L/4,(C1)+R , 0, m};
Point(27) = {3*L/4,(C1)-R , 0, m};

Circle(25) = {26, 25, 27};
Circle(26) = {27, 25, 26};

Point(28) = {3*L/4,2*C1 , 0, m};
Point(29) = {3*L/4,2*(C1)+R , 0, m};
Point(210) = {3*L/4,2*(C1)-R , 0, m};

Circle(27) = {210, 28, 29};
Circle(28) = {29, 28, 210};

Point(211) = {3*L/4,3*C1 , 0, m};
Point(212) = {3*L/4,3*(C1)+R , 0, m};
Point(213) = {3*L/4,3*(C1)-R , 0, m};
 
Circle(29) = {213, 211, 212};
Circle(210) = {212, 211, 213};


Point(35) = {1.5*L/4,1.5*C1 , 0, m};
Point(36) = {1.5*L/4,1.5*(C1)+R , 0, m};
Point(37) = {1.5*L/4,1.5*(C1)-R , 0, m};

Circle(35) = {36, 35, 37};
Circle(36) = {37, 35, 36};
 


Point(45) = {1.5*L/4,2.5*C1 , 0, m};
Point(46) = {1.5*L/4,2.5*(C1)+R , 0, m};
Point(47) = {1.5*L/4,2.5*(C1)-R , 0, m};

Circle(45) = {46, 45, 47};
Circle(46) = {47, 45, 46};

Point(55) = {2.5*L/4,1.5*C1 , 0, m};
Point(56) = {2.5*L/4,1.5*(C1)+R , 0, m};
Point(57) = {2.5*L/4,1.5*(C1)-R , 0, m};

Circle(55) = {56, 55, 57};
Circle(56) = {57, 55, 56};
 


Point(65) = {2.5*L/4,2.5*C1 , 0, m};
Point(66) = {2.5*L/4,2.5*(C1)+R , 0, m};
Point(67) = {2.5*L/4,2.5*(C1)-R , 0, m};

Circle(65) = {66, 65, 67};
Circle(66) = {67, 65, 66};
//+
Line Loop(211) = {3, 4, 1, 2};
//+
Line Loop(212) = {9, 10};
//+
Line Loop(213) = {7, 8};
//+
Line Loop(214) = {6, 5};
//+
Line Loop(215) = {45, 46};
//+
Line Loop(216) = {36, 35};
//+
Line Loop(217) = {110, 19};
//+
Line Loop(218) = {18, 17};
//+
Line Loop(219) = {16, 15};
//+
Line Loop(220) = {65, 66};
//+
Line Loop(221) = {55, 56};
//+
Line Loop(222) = {210, 29};
//+
Line Loop(223) = {28, 27};
//+
Line Loop(224) = {25, 26};
//+
Plane Surface(225) = {211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224};
//+
Transfinite Line {3, 1, 4, 2, 9, 10, 8, 7, 6, 5, 46, 45, 36, 35, 110, 19, 17, 18, 16, 15, 65, 66, 56, 55, 28, 27, 210, 29, 25, 26} = 100 Using Progression 1;
//+
Transfinite Surface {225};
//+
Recombine Surface {225};


Physical Line(0) = {4, 1, 2, 3, 9, 10, 8, 7, 6, 5, 45, 46, 36, 35, 15, 16, 18, 17, 110, 19, 65, 66, 56, 55, 28, 27, 210, 29, 25, 26};
//+
Physical Surface(227) = {225};
