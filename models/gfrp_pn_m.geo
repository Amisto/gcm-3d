cl1 = 1;
Point(1) = {-4.50, 0.1080, -4.50, cl1};
Point(2) = {-4.50, 0.3220, -4.50, cl1};
Point(3) = {4.50, 0.3220, -4.50, cl1};
Point(4) = {4.50, 0.1080, -4.50, cl1};
Point(5) = {4.50, 0.1080, 4.5, cl1};
Point(6) = {-4.50, 0.1080, 4.5, cl1};
Point(7) = {-4.50, 0.3220, 4.5, cl1};
Point(8) = {4.50, 0.3220, 4.5, cl1};
Line(1) = {2, 3};
Line(2) = {3, 8};
Line(3) = {8, 7};
Line(4) = {7, 2};
Line(5) = {2, 1};
Line(6) = {1, 4};
Line(7) = {4, 5};
Line(8) = {5, 6};
Line(9) = {6, 1};
Line(10) = {6, 7};
Line(11) = {5, 8};
Line(12) = {4, 3};
Line Loop(14) = {6, 12, -1, 5};
Plane Surface(14) = {14};
Line Loop(16) = {4, 5, -9, 10};
Plane Surface(16) = {16};
Line Loop(18) = {8, 10, -3, -11};
Plane Surface(18) = {18};
Line Loop(20) = {11, -2, -12, 7};
Plane Surface(20) = {20};
Line Loop(22) = {3, 4, 1, 2};
Plane Surface(22) = {22};
Line Loop(24) = {9, 6, 7, 8};
Plane Surface(24) = {24};
Surface Loop(26) = {22, 18, 24, 16, 14, 20};
Volume(26) = {26};
