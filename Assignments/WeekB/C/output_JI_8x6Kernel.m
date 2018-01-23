your_version = 'JI\_8x6Kernel';
 
% number of repeats:% 3
% enter first, last, inc:% 48 960 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   960 4.8000e-02 3.6864e+01    1.5316e-01 1.1553e+01 9.6634e-13
   912 4.1442e-02 3.6608e+01    1.0075e-01 1.5057e+01 1.1084e-12
   864 3.5166e-02 3.6682e+01    9.8383e-02 1.3111e+01 8.8107e-13
   816 2.9646e-02 3.6655e+01    7.2724e-02 1.4943e+01 8.2423e-13
   768 2.4342e-02 3.7218e+01    7.4910e-02 1.2094e+01 7.3896e-13
   720 2.0016e-02 3.7294e+01    5.1195e-02 1.4582e+01 7.1054e-13
   672 1.6382e-02 3.7048e+01    4.5252e-02 1.3412e+01 6.5370e-13
   624 1.3234e-02 3.6720e+01    3.4052e-02 1.4270e+01 5.4001e-13
   576 1.0846e-02 3.5240e+01    2.9151e-02 1.3111e+01 5.1159e-13
   528 8.4041e-03 3.5030e+01    1.7303e-02 1.7014e+01 4.5475e-13
   480 6.1855e-03 3.5758e+01    1.2906e-02 1.7138e+01 3.2685e-13
   432 4.5587e-03 3.5370e+01    9.4949e-03 1.6982e+01 2.8422e-13
   384 3.3061e-03 3.4253e+01    8.6029e-03 1.3164e+01 2.1316e-13
   336 2.2320e-03 3.3990e+01    4.4705e-03 1.6970e+01 1.7053e-13
   288 1.4513e-03 3.2919e+01    2.7992e-03 1.7068e+01 1.1369e-13
   240 8.1924e-04 3.3748e+01    1.6345e-03 1.6915e+01 4.2633e-14
   192 4.3139e-04 3.2814e+01    8.0652e-04 1.7552e+01 2.8422e-14
   144 1.9810e-04 3.0147e+01    2.0165e-04 2.9615e+01 2.8422e-14
    96 6.8550e-05 2.5813e+01    6.0178e-05 2.9404e+01 1.4211e-14
    48 1.4080e-05 1.5709e+01    7.0290e-06 3.1467e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 1.108447e-12.
