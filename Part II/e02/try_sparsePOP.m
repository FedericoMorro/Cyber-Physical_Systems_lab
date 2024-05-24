clear
close all
clc

format compact

objPoly.typeCone = 1;
objPoly.dimVar = 3;
objPoly.degree = 1;
objPoly.noTerms = 3;
objPoly.supports = [
    1   0   0
    0   1   0
    0   0   1
];
objPoly.coef = [
    -2
    3
    -2
];


ineqPolySys{1}.typeCone = 1;
ineqPolySys{1}.dimVar = 3;
ineqPolySys{1}.degree = 2;
ineqPolySys{1}.noTerms = 8;
ineqPolySys{1}.supports = [
    0   0   0
    1   0   0
    0   1   0
    0   0   1
    2   0   0
    0   2   0
    0   1   1
    0   0   2
];
ineqPolySys{1}.coef = [
    19
    -17
    8
    -14
    6
    3
    -2
    3
];

ineqPolySys{2}.typeCone = 1;
ineqPolySys{2}.dimVar = 3;
ineqPolySys{2}.degree = 1;
ineqPolySys{2}.noTerms = 4;
ineqPolySys{2}.supports = [
    0   0   0
    1   0   0
    0   1   0
    0   0   1
];
ineqPolySys{2}.coef = [
    5
    -1
    -2
    -1
];

ineqPolySys{3}.typeCone = -1;
ineqPolySys{3}.dimVar = 3;
ineqPolySys{3}.degree = 1;
ineqPolySys{3}.noTerms = 3;
ineqPolySys{3}.supports = [
    0   0   0
    0   1   0
    0   0   1
];
ineqPolySys{3}.coef = [
    -7
    5
    2
];

ubd = [2 1 1e10];
lbd = [0 0 -1e10];

param.relaxOrder = 1;
param.POPsolver = 'active-set';

[param,SDPobjValue,POP,elapsedTime,SDPsolverInfo,SDPinfo] = ...
    sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
POP.objValueL
POP.xVectL