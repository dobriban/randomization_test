function T = t2(X1,X2,n1,n2)

Del = (mean(X1,1)-mean(X2,1));
S2 = ((n1-1)*var(X1)+ (n2-1)*var(X2))/(n1+n2-2);
T = Del/sqrt(S2*(1/n1+1/n2));