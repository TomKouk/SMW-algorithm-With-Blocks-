A = MxMake_1781('mc',1024);
for k=1:1024/2
    x(2*k-1)=1;
    x(2*k)=(-1)^(k+1)*1/(2*k);
end
x=x';
b=A*x;

t1=0;
for j=1:4
  x=SMW_solve_1781(A,b,4,5,6,'colwise');
  tic
  x=SMW_solve_1781(A,b,4,5,6,'colwise');
  toc;
  t1=t1+toc;
end
t1=t1/4;

tries=[1,2,4,8,16];
for i=1:5
  t2(i)=0;
  x=SMW_solve_1781_block(A,b,4,5,6,tries(i),'colwise');
  for j=1:4
    tic
    x=SMW_solve_1781_block(A,b,4,5,6,tries(i),'colwise');
    t2(i)=t2(i)+toc;
  end
  t2(i)=t2(i)/4;
end
