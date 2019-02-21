function T = T_matrix(a)

T=[cosd(a) sind(a) 0 0 0 0
    -sind(a) cosd(a) 0 0 0 0
    0 0 1 0 0 0
    0 0 0 cosd(a) sind(a) 0
    0 0 0 -sind(a) cosd(a) 0
    0 0 0 0 0 1];

    end