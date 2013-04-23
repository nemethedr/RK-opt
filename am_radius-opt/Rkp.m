function [R,alpha,beta]=Rkp(k,p)
% function [R,alpha,beta]=Rkp(k,p)
%
% Find the optimal SSP k-step explicit LMM with order of accuracy p.
%
% Inputs: 
%       * `k` = # of steps
%       * `p` = order of accuracy
%
% Outputs: 
%       * `\alpha, \beta` = the coefficients of the method
%
% Requires MATLAB's optimization toolbox for the LP solver.
    
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

M=2*k;              %Number of decision variables (coefficients)
rmax=2.0001; rmin=0;  %Upper and lower bounds for R
r=rmax;               %Initial guess
c=zeros(M,1); d=zeros(p+1,1); B=zeros(p+1,M);
%=========================================================
  
%=========================================================
%Find R by bisection
while (rmax-rmin>acc) 
  %Set up equality constraints
  %g: First k+1 unkowns are beta's, last k are gammas
  for i=0:p
    d(i+1)=1;
    for j=0:k-1
      if (i+j==0) %Avoid divide-by-zero
        B(i+1,j+1)=r/k^i;
      else
        B(i+1,j+1)=(r*j + i)*(j/k)^(i-1) / k;
      end %if
      %B(i+1,k+j+2)=(j/k)^i;
      B(i+1,k+j+1)=(j/k)^i;
    end
    %B(i+1,k+1)=i/k;
  end
  %Test feasibility for this value of r
  [x,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts);
  if exitflag==1
    rmin=r; r=(r+rmax)/2;
  else
    rmax=r; r=(rmin+r)/2;
  end
end

%=========================================================
%Now get a feasible solution so we have the coefficients
r=rmin;
for i=0:p
  d(i+1)=1;
  for j=0:k-1
    if (i+j==0) %Avoid divide-by-zero
      B(i+1,j+1)=r/k^i;
    else
      B(i+1,j+1)=(r*j + i)*(j/k)^(i-1) / k;
    end %if
    %B(i+1,k+j+2)=(j/k)^i;
    B(i+1,k+j+1)=(j/k)^i;
  end
  %B(i+1,k+1)=i/k;
end
[g,lambda,exitflag]=linprog(c,[],[],B,d,zeros(M,1),zeros(M,1)+1.e6,c,opts);
%=========================================================
% Solve the dual problem
% A*y <= b y_i: i=0...p, j=0...(k-1),0...(k-1),k
% epsilon=minimum value of the polynomial at k
% average=the average value of the polynomial
r=rmax;
epsilon=0.0; %not used
average=1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Aeq=zeros(1,p+1);
beq=average;
A=zeros(2*k+1,p+1);
b=zeros(2*k+1,1); b(end)=-epsilon;
y=zeros(1,p+1);
f=zeros(1,p+1);
for j=0:k-1
  for i=0:p
    A(j+1,i+1)=-j^i;
    A(j+k+1,i+1)=-r*j^i;
    if (i > 0)
      A(j+k+1,i+1)=A(j+k+1,i+1)-i*j^(i-1);
    end
  end
end
for i=0:p
  A(2*k+1,i+1)=k^i;
end
for i=0:p
  Aeq(1,i+1)=-sum(A(1:k,i+1))/k;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[y,lambda,exitflag]=linprog(f,A,b,Aeq,beq,[ ],[ ],y,opts);
ypoly=y(end:-1:1);% matlab ordering
rts=roots(ypoly);

%=========================================================
% Prepare outputs
R=rmin;
beta=g(1:k); alpha=g(k+1:end)+R*beta;
g
d;
B(:,[1 2 6]);
B;
if (rmin > 0.0)
  y
  rts
end
