//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-0.7, -0.3, 0, 1, 1, 0};
MeshSize {1, 4, 3, 2} = 1;
//+
//Point(5) = {-0., 0., 0, 1.0};
//Point{5} In Surface {1};
//+
Physical Curve("clamped", 5) = {4};
//+
Physical Surface("bulke", 6) = {1};
//+
//Physical Point("forcex", 7) = {3};
//+
//Physical Line("forcey", 8) = {2};
//+
