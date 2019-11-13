cases={'had','trihad','toep','mc','wathen','CollegeMsg'};
n=[64,64,64,400,12,1];

i=0;
for i=1:6
        
    A=MxMake_1781(cases{i},n(i));
    cond_num(i)=cond(A,inf);
    
    if i==5, n(i)=443; end
    x=zeros(n(i),1);
    for k=1:n(i)/2
        x(2*k-1)=1;
        x(2*k)=(-1)^(k+1)*1/(2*k);
    end    
    b=A*x;
        
    x1=SMW_solve_1781(A,b,1,2,3,'colwise');
    back_error1(i)=normest(A*x1-b)/(normest(A)*normest(x1)+normest(b));
    forwar_error1(i)=2*back_error1(i)*cond_num(i);
    actual_error1(i)=normest(x1-x)/normest(x);
    
    x2=A\b; 
    back_error2(i)=normest(A*x2-b)/(normest(A)*normest(x2)+normest(b));
    forwar_error2(i)=2*back_error2(i)*cond_num(i);
    actual_error2(i)=normest(x-x2)/normest(x);
    
    x=[];
end
