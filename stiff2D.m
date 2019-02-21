function Frame2D = stiff2D(E,I,A,L)

a=(E*A)/L;
b=(E*I)/L;
c=(E*I)/L^2;
d=(E*I)/L^3;

Frame2D=[a 0 0 -a 0 0
         0 12*d 6*c 0 -12*d 6*c
         0 6*c 4*b 0 -6*c 2*b
         -a 0 0 a 0 0
         0 -12*d -6*c 0 12*d -6*c
         0 6*c 2*b 0 -6*c 4*b];

end