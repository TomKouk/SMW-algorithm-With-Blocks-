function x = SMW_solve_1781(A,b,M,P,Q,sdir)

[~,n]=size(A);
A0=diag(diag(A));

%colwise case
if strcmp(sdir,'colwise'), 
  P=A0\(A-A0);
  Q=eye(n);
  
%rowwise case  
elseif strcmp(sdir,'rowwise'), 
  P=eye(n);
  Q=(A0\(A-A0))';
  
%other cases  
else
  %BLAS-3
  x=M^(-1)*b-M^(-1)*(P*((I+Q'*M^(-1)*P)^(-1)*(Q'*(M^(-1)*b))));
  return
end

x0=A0\b;
for i=1:(n-1)
  dot1=Q(:,i)'*x0; %DOT
  dot2=Q(:,i)'*P(:,i); %DOT
  num1=dot1/(1+dot2); %NUMs
  x0=x0-num1*P(:,i); %SAXPY
  for j=i+1:n
    dot3=Q(:,i)'*P(:,j); %DOT
    dot4=Q(:,i)'*P(:,i); %DOT
    num2=dot3/(1+dot4); %NUMs
    P(:,j)=P(:,j)-P(:,i)*num2; %SAXPY
  end
end
dot5=Q(:,end)'*x0; %DOT
dot6=Q(:,end)'*P(:,end); %DOT
num3=dot5/(1+dot6); %NUMs
x=x0-num3*P(:,end); %SAXPY

end

