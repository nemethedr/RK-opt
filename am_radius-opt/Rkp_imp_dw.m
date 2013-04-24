function [R,alpha,beta,tbeta]=Rkp_imp_dw(k,p)
% function [R,alpha,beta]=Rkp_imp_dw(k,p)
%
% Finds the optimal k-step implicit LMM with order of accuracy p
% allowing downwinding
%
% Inputs: k = # of steps
%        p = order of accuracy
%
% Outputs: alpha, beta, tbeta = the coefficients of the method
%    
% Depends on MATLAB's optimization toolbox for the LP solver

%=========================================================
%Initialize
clear tbtest;
tbtest=which('linprog');
if ~exist(tbtest)
  disp('MATLAB optimization toolbox is required')
  return
end
%Set options for linprog
opts=optimset('TolX',1.e-15,'TolFun',1.e-15,'MaxIter',10000000,...
               'LargeScale','off','Simplex','off','Display','off');
acc=1.e-15; %Accuracy of bisection search

M=3*k+2;              %Number of decision variables (coefficients)
rmax=M+.0001; rmin=0; %Upper and lower bounds for R
r=rmax;               %Initial guess
c=zeros(M,1); d=zeros(p+1,1); B=zeros(p+1,M);
lb=zeros(M,1); ub=zeros(M,1)+1.e6;
%=========================================================

 
while (rmax-rmin>acc) %Find R by bisection
 %Set up equality constraints
 %g: First k+1 unkowns are beta's, next k+1 are \tilde{\beta}'s, last k are gammas
  for i=0:p
    d(i+1)=k^i;
    for j=0:k-1
      if (i+j==0) %Avoid divide-by-zero
        B(i+1,  j+1)=r;                       %beta
        B(i+1,k+1+j+1)=r;                    %tbeta
      else
        B(i+1,  j+1)=r*j^i + i*j^(i-1);       %beta
        B(i+1,k+1+j+1)=r*j^i - i*j^(i-1);    %tbeta
      end %if
      B(i+1,2*(k+1)+j+1)=j^i;                 %gamma
    end
    B(i+1,k+1)=i*k^(i-1);
    B(i+1,2*(k+1))=i*k^(i-1);
  end
  %Test feasibility for this value of r
  [g,lambda,exitflag]=linprog(c,[],[],B,d,lb,ub,c,opts);
  if exitflag==1
    rmin=r; r=(r+rmax)/2;
  else
    rmax=r; r=(rmin+r)/2;
  end
end

%Now get a feasible solution so we have the coefficients
r=rmin;
for i=0:p
  d(i+1)=k^i;
  for j=0:k-1
    if (i+j==0) %Avoid divide-by-zero
      B(i+1,j+1)=r;
      B(i+1,k+1+j+1)=r;
    else
      B(i+1,j+1)=r*j^i + i*j^(i-1);
      B(i+1,k+1+j+1)=r*j^i - i*j^(i-1);
    end %if
    B(i+1,2*(k+1)+j+1)=j^i;
  end
  B(i+1,k+1)=i*k^(i-1);
  B(i+1,2*(k+1))=i*k^(i-1);
end
[g,lambda,exitflag]=linprog(c,[],[],B,d,lb,ub,c,opts);

%=========================================================
% Solve the dual problem
% A*y <= b y_i: i=0...p, j=0...(k-1),0...(k-1),0...(k-1),k
% epsilon=minimum value of the polynomial at k
% average=the average value of the polynomial
r=rmax;
epsilon=0.0; %not used
average=1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Aeq=zeros(2,p+1);
beq=zeros(2,1); beq(end)= average;
A=zeros(3*k+3,p+1);
b=zeros(3*k+3,1); b(end)=-epsilon;
y=zeros(1,p+1);
f=zeros(1,p+1);
for j=0:k-1
  for i=0:p
    A(j+1,i+1)=-j^i;
    A(j+k+1,i+1)=-r*j^i;
    A(j+2*k+1,i+1)=-r*j^i;
    if (i > 0)
      A(j+k+1,i+1)=A(j+k+1,i+1)-i*j^(i-1);
      A(j+2*k+1,i+1)=A(j+2*k+1,i+1)+i*j^(i-1);
    end
  end
end
for i=0:p
  A(3*k+1,i+1)=k^i;
end
for i=0:p
  Aeq(1,i+1)=i*k^(i-1);
  Aeq(2,i+1)=-sum(A(1:k,i+1))/k;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[y,lambda,exitflag]=linprog(f,A,b,Aeq,beq,[ ],[ ],y,opts);
ypoly=y(end:-1:1);% matlab ordering
rts=roots(ypoly);
%=========================================================
R=rmin;
beta=g(1:k+1); tbeta=g(k+2:2*k+2); alpha=g(2*k+3:end)+R*(beta(1:k)+tbeta(1:k));
%alpha=wrev(alpha); beta=wrev(beta); tbeta=wrev(tbeta);
alpha=alpha(end:-1:1); beta=beta(end:-1:1); tbeta=tbeta(end:-1:1);
if (R > 0.0)
  y
  rts
end

