function A = MxMake_1781(mx_id,n)

if strcmp(mx_id,'had'),A=hadamard(n);end
if strcmp(mx_id,'trihad'),A=triu(hadamard(n));end
if strcmp(mx_id,'toep'),A=toeplitz([4,-1,zeros(1,n-2)]);end
if strcmp(mx_id,'mc'),
    for i=1:n
        for j=1:n
            if i==j,A(i,i)=1+i;else
               A(i,j)=1/abs(i+j)^2;end
        end
     end
end

if strcmp(mx_id,'wathen'),A=gallery('wathen',n,n-1);end

if strcmp(mx_id,'CollegeMsg'),
  websave('collegeMsg','https://sparse.tamu.edu/mat/SNAP/CollegeMsg.mat');
  load('CollegeMsg.mat');
  A=eye(1899)-0.85*Problem.A;
end

end

