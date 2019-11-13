function x = SMW_solve_1781_block(A,b,M,P,Q,blk,sdir)

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
for i=1:blk:(n-blk)
  mm1=Q(:,i:i+blk-1)'*x0; %matrix*vector
  mm2=Q(:,i:i+blk-1)'*P(:,i:i+blk-1); %matrix*matrix
  inver=(eye(blk)+mm2)^(-1); %inversion

  x0=x0-P(:,i:i+blk-1)*(inver*mm1);
    for step1=(i+blk):blk:n
      step2=step1+blk-1;
      if (step1+blk-1)>n, step2=n; end

      mm1=Q(:,i:i+blk-1)'*P(:,i:i+blk-1); %matrix*matrix
      mm2=Q(:,i:i+blk-1)'*P(:,step1:step2); %matrix*matrix
      inver=(eye(blk)+mm1)^(-1); %inversion

      P(:,step1:step2)=P(:,step1:step2)-P(:,i:i+blk-1)*inver*mm2;
    end
end
step3=n-blk+1;
if mod(n,blk)>0, step3=n-mod(n,blk)+1; end

mm1=Q(:,step3:n)'*P(:,step3:n); %matrix*matrix
inver=(eye(n-step3+1)+mm1)^(-1); %inversion
mm2=Q(:,step3:n)'*x0; %matrix*vector

x=x0-P(:,step3:n)*inver*mm2;

end

