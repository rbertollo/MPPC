function [A,C,G] = getMatricesEKF(x)
lambda = 0.1;


A = [          -x(2)*x(7),                -x(1)*x(7),                      0,              0,     0,     0,      -x(2)*x(1),   0,     0,   0,    0;
               x(2)*x(7),       x(1)*x(7) - lambda - x(8),                0,              0,     0,     0,      x(2)*x(1),  -x(2),   0,   0,    0;
                       0,                    x(8),                  -x(9) - lambda,       0,     0,     0,          0,       x(2), -x(3), 0,    0;
                       0,                     0,                     x(9),       -x(10) - x(11), 0,     0,          0,        0,    x(3), -x(4), -x(4);
                       0,                   lambda,                lambda,              x(10),   0,     0,          0,        0,     0,  x(4),  0;
                       0,                      0,                      0,                x(11),  0,     0,          0,        0,     0,    0, x(4);
                       0,                      0,                      0,                 0,     0,     0,          0,        0,    0,    0,    0;
                       0,                      0,                      0,                 0,     0,     0,          0,        0,    0,    0,    0;
                       0,                      0,                      0,                 0,     0,     0,          0,        0,    0,    0,    0;
                       0,                      0,                      0,                 0,     0,     0,          0,        0,    0,    0,    0;
                       0,                      0,                      0,                 0,     0,     0,          0,        0,    0,    0,    0;
];


C = [ diag([1 1 1 1 1 1]), zeros(6,5)];
G = [ eye(length(x))];
end