function [R,alpha,beta,tbeta]=Rkp_dw(k,p)
% function [R,alpha,beta,tbeta]=Rkp_dw(k,p)
%
% Finds the optimal SSP k-step explicit LMM with order of accuracy p
% allowing downwind operators
%
% Inputs: k = # of steps
%         p = order of accuracy
%
% Outputs: alpha, beta, tbeta = the coefficients of the method
%
% The method is given by
% `u_n = \sum_{j=0}^{k-1} (\alpha[j] + \beta[j] F(u_{n-k+j} + tbeta[j] tF(u_{n-k+j}))`
% where tF(u) is the negated downwind operator.
%
% Depends on MATLAB's optimization toolbox for the LP solver
    
%=========================================================
%Set options for linprog
opts=optimset('TolX',1.e-15,'TolFun',1.e-15,'MaxIter',10000000,...
               'LargeScale','off','Simplex','off','Display','off');
acc=1.e-15; %Accuracy of bisection search
%=========================================================

clear B d;
M=3*k;                %Number of decision variables (coefficients)
rmax=3*k; rmin=0;     %Upper and lower bounds for R
r=rmax;               %Initial guess
c=zeros(M,1);
lb=zeros(M,1); ub=zeros(M,1)+1.e6; %Lower and upper bounds on coefficients

%=========================================================
%Find R by bisection
while (rmax-rmin>acc) 
  %Set up equality constraints
  %g: First k unkowns are beta's, next k are \tilde{\beta}'s, last k are gammas
  for i=0:p
    d(i+1)=1;
    for j=0:k-1
      if (i+j==0) %Avoid divide-by-zero
        B(i+1,  j+1)=r;
        B(i+1,k+j+1)=r;
      else
        B(i+1,  j+1)=r*(j/k)^i + i/k*(j/k)^(i-1);
        B(i+1,k+j+1)=r*(j/k)^i - i/k*(j/k)^(i-1);
      end %if
      B(i+1,2*k+j+1)=(j/k)^i;
    end
  end
  %Test feasibility for this value of r
  [g,lambda,exitflag]=linprog(c,[],[],B,d,lb,ub,c,opts);
  if exitflag==1
    rmin=r; r=(r+rmax)/2;
  else
    rmax=r; r=(rmin+r)/2;
  end
end
%=========================================================

%=========================================================
%Now get a feasible solution so we have the coefficients
r=rmin;
for i=0:p
  d(i+1)=1;
  for j=0:k-1
    if (i+j==0) %Avoid divide-by-zero
      B(i+1,j+1)=r;
      B(i+1,k+j+1)=r;
    else
      B(i+1,  j+1)=r*(j/k)^i + i/k*(j/k)^(i-1);
      B(i+1,k+j+1)=r*(j/k)^i - i/k*(j/k)^(i-1);
    end %if
    B(i+1,2*k+j+1)=(j/k)^i;
  end
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
Aeq=zeros(1,p+1);
beq=average;
A=zeros(3*k+1,p+1);
b=zeros(3*k+1,1); b(end)=-epsilon;
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
  Aeq(1,i+1)=-sum(A(1:k,i+1))/k;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[y,lambda,exitflag]=linprog(f,A,b,Aeq,beq,[ ],[ ],y,opts);
ypoly=y(end:-1:1);% matlab ordering
rts=roots(ypoly);
%=========================================================

%Prepare outputs
R=rmin;
beta=g(1:k); tbeta=g(k+1:2*k); alpha=g(2*k+1:end)+R*(beta+tbeta);
alpha=alpha(end:-1:1); beta=beta(end:-1:1); tbeta=tbeta(end:-1:1);
if (R > 0.0)
  y
  rts
end
