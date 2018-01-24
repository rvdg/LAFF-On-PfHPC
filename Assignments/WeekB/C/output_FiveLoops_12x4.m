your_version = 'GemmFiveLoops\_12x4';
 
% number of repeats:% 3
% enter first, last, inc:% 48 960 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   960 4.7924e-02 3.6923e+01    7.3140e-02 2.4193e+01 9.6634e-13
   912 4.1315e-02 3.6720e+01    6.2280e-02 2.4359e+01 1.1084e-12
   864 3.5091e-02 3.6760e+01    5.3600e-02 2.4066e+01 8.8107e-13
   816 2.9668e-02 3.6627e+01    4.6984e-02 2.3129e+01 8.2423e-13
   768 2.4394e-02 3.7139e+01    4.4995e-02 2.0135e+01 7.3896e-13
   720 2.0052e-02 3.7228e+01    3.0606e-02 2.4390e+01 7.1054e-13
   672 1.6385e-02 3.7042e+01    2.4850e-02 2.4424e+01 6.5370e-13
   624 1.3224e-02 3.6747e+01    2.0296e-02 2.3943e+01 5.4001e-13
   576 1.0490e-02 3.6435e+01    1.6692e-02 2.2898e+01 5.1159e-13
   528 8.1648e-03 3.6057e+01    1.2075e-02 2.4381e+01 4.5475e-13
   480 6.0291e-03 3.6686e+01    9.0948e-03 2.4320e+01 3.2685e-13
   432 4.4489e-03 3.6243e+01    6.6582e-03 2.4217e+01 2.8422e-13
   384 3.2129e-03 3.5247e+01    5.5083e-03 2.0559e+01 2.1316e-13
   336 2.1785e-03 3.4826e+01    3.0375e-03 2.4976e+01 1.7053e-13
   288 1.4185e-03 3.3680e+01    1.8345e-03 2.6043e+01 1.1369e-13
   240 8.0541e-04 3.4328e+01    1.0235e-03 2.7014e+01 4.2633e-14
   192 4.2677e-04 3.3169e+01    5.1852e-04 2.7301e+01 2.8422e-14
   144 1.9634e-04 3.0417e+01    2.1615e-04 2.7629e+01 2.8422e-14
    96 6.7528e-05 2.6204e+01    6.2377e-05 2.8367e+01 1.4211e-14
    48 1.3900e-05 1.5913e+01    7.3540e-06 3.0077e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 1.108447e-12.
