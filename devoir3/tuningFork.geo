

lc = 0.01;
a = 1e-2;

Point(3) = {3.81*a, 0, 0, lc};
Point(4) = {21.9075*a, 0, 0, lc};
Point(5) = {21.9075*a, 1.27*a, 0, lc};
Point(6) = {6.35*a, 1.27*a, 0, lc};
Point(7) = {6.35*a, 2.54*a, 0, lc};
Point(8) = {21.9075*a, 2.54*a, 0, lc};
Point(9) = {21.9075*a, 3.81*a, 0, lc};
Point(10) = {3.81*a, 3.81*a, 0, lc};

Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,3};

Curve Loop(1) = {3, 4, 5, 6, 7, 8, 9, 10};
Plane Surface(1) = {1};
Physical Surface("bulke", 7) = {1};
Physical Curve("clamped", 5) = {10};
