cl1=0.03;
Point(1) = {0.5,0,0,cl1};
Point(2) = {0,0.5,0,cl1};
Point(3) = {-0.5,0,0,cl1};
Point(4) = {0,-0.5,0,cl1};
Point(5) = {0,0,0.5,cl1};
Point(6) = {0,0,-0.5,cl1};
Point(7) = {0,0,0,cl1};

Circle(1) = {3, 7, 6};
Circle(2) = {2, 7, 6};
Circle(3) = {1, 7, 6};
Circle(4) = {4, 7, 6};
Circle(5) = {3, 7, 5};
Circle(6) = {3, 7, 4};
Circle(7) = {4, 7, 1};
Circle(8) = {1, 7, 2};
Circle(9) = {2, 7, 3};
Circle(10) = {1, 7, 5};
Circle(11) = {2, 7, 5};
Circle(12) = {5, 7, 4};

Line Loop(13) = {11, -10, 8};
Ruled Surface(14) = {13};


Line Loop(15) = {11, -5, -9};
Ruled Surface(16) = {-15};


Line Loop(17) = {9, 1, -2};
Ruled Surface(18) = {-17};


Line Loop(19) = {2, -3, 8};
Ruled Surface(20) = {-19};


Line Loop(21) = {4, -3, -7};
Ruled Surface(22) = {21};


Line Loop(23) = {1, -4, -6};
Ruled Surface(24) = {23};


Line Loop(25) = {7, 10, 12};
Ruled Surface(26) = {25};


Line Loop(27) = {12, -6, 5};
Ruled Surface(28) = {-27};

Dilate {{0, 0, 0}, 0.9} {
  Duplicata { Surface{18, 20, 14, 16, 28, 24, 22, 26}; }
}
Dilate {{0, 0, 0}, 0.8} {
  Duplicata { Surface{18, 20, 14, 16, 28, 24, 22, 26}; }
}
Dilate {{0, 0, 0}, 0.7} {
  Duplicata { Surface{18, 20, 14, 16, 28, 24, 22, 26}; }
}
Dilate {{0, 0, 0}, 0.6} {
  Duplicata { Surface{18, 20, 14, 16, 28, 24, 22, 26}; }
}

Physical Surface(1) = {14, 16, 18, 20, 22, 24, 28, 26};
Physical Surface(2) = {37, 41, 29, 33, 53, 49, 45, 57};
Physical Surface(3) = {66, 70, 58, 62, 82, 78, 74, 86};
Physical Surface(4) = {95, 99, 87, 91, 111, 107, 103, 115};
Physical Surface(5) = {124, 128, 116, 120, 140, 136, 132, 144};

Surface Loop(1) = {14, 16, 18, 20, 22, 24, 28, 26};
Surface Loop(2) = {37, 41, 29, 33, 53, 49, 45, 57};
Surface Loop(3) = {66, 70, 58, 62, 82, 78, 74, 86};
Surface Loop(4) = {95, 99, 87, 91, 111, 107, 103, 115};
Surface Loop(5) = {124, 128, 116, 120, 140, 136, 132, 144};

Volume(1) = {1, 2};
Volume(2) = {2, 3};
Volume(3) = {3, 4};
Volume(4) = {4, 5};
Volume(5) = {5};

Physical Volume(1) = {1};
Physical Volume(2) = {2};
Physical Volume(3) = {3};
Physical Volume(4) = {4};
Physical Volume(5) = {5};
