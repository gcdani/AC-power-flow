% Type:
% 1 - PQ
% 2 - PV
% 3 - slack

%         Bus      Type V    Degree Pg     Qg     Pd     Qd   Gs   Bs  Qmin    Qmax
bus  = [   1        3   1.0   0     0      0      0       0     0    0    0       0;
           2        2   1.0   0   150.0   100    200.0   100.0  0    0   -50      100];
       

%        from to   R       X     b      a    shift
line = [  1   2    0.020   0.20  0.50    1      0;
          1   2    0.025   0.25  0.45    1      0];