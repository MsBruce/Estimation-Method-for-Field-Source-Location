
clear;clc
load data.mat
figure
contourf(T90)
%%
[ BNSS ] = New_NSS(T90,201,201,-4000,4000,-4000,4000,0);
figure
contourf(BNSS)
title('BNSS')


