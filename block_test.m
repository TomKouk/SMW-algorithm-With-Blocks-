A = MxMake_1781('mc',64);
for k=1:64/2
    x(2*k-1)=1;
    x(2*k)=(-1)^(k+1)*1/(2*k);
end
x=x';
b=A*x;
k=cond(A,inf); 
tries=[2,5,15,24,68];
for i=1:5
  block_x = SMW_solve_1781_block(A,b,4,5,6,tries(i),'colwise');
	rigal1(i)=normest(A*block_x-b)/(normest(A)*normest(block_x)+normest(b));
	rigal2(i)=2*rigal1(i)*k;
end


