function Steepest_Descent(x0, epsilon, A0, b0, n)
% Steepest_Descent.m, (c) Akash Varma, 2021
%
% created: 	Tue May 05 2021 
% author:  	Akash Varma
% email:   	a1678314@student.adelaide.edu.au
% 
% Steepest_Descent Performs the steepest descent algorithm to optimise a
% quadratic defined by two matrices A and b.
%     A,b are the matrices defining the function to be optimised; 
%     x(i,:) is the value of xk in the ith iteration; 
%     u(k,:) is the direction used in the kth iteration; 
%     lam(k) gives the value of the step size in the kth iteration; 
%     convk(k,:) gives the convergence rate and its upper bound in the 
%     kth iteration; k is the iteration counter;
%
%
% INPUTS:
%        x0 = The starting position from which the algorithm commences;
%     A0,b0 = The nxn and nx1 matrices that define the quadratic to be optimised; 
%     epsilon = The upper bound on the magnitude of the gradient of the final estimate.
%     n = The dimensions of matrix A;


%Initialise inputs
A = A0;
b = b0;

%Compute the condition number of A
s = svd(A);
r = s(1)/s(n);

%convergence bounds
delta = ((r-1)/(r+1));

%Compute the true minimiser of q(x)
true_min = (-A\b)';

%initialise iteration counter
k = 1;

%setup matrices for storing vectors values for each iteration
x = zeros(100, n);
u = zeros(100, n);

%setup arrays for storing scalars for each iteration
lam = zeros(100, 1);

%compute and store the zeroth values of x,g,u
x(1,:) = x0';
u(1,:) = -(A*x(1,:)' + b)';

%magic number for entering loop
norm_uk = 100;

%perform the algorithm 
while norm_uk >= epsilon
    
    u(k,:) = -(A*x(k,:)' + b)';
    
    lam(k) = (u(k,:)*u(k,:)')/(u(k,:)*A*u(k,:)');

    x(k+1,:) = x(k,:) + lam(k)*u(k,:);
    
    norm_uk = norm(u(k,:));
    
    k = k + 1;
end

%Convergence of each iteration
convk = zeros(k-1,2);

%Gap between q(x0) and q(x*)
initgap = ((1/2)*x(1,:)*A*x(1,:)' + b'*x(1,:)') - ((1/2)*true_min*A*true_min' + b'*true_min');

%Compute convergence rate for each iteration
for i=1:1:k
    
    gapi = ((1/2)*x(i,:)*A*x(i,:)' + b'*x(i,:)') - ((1/2)*true_min*A*true_min' + b'*true_min');
    
    convk(i,1) = gapi;
    convk(i,2) = delta^(2*(i-1))*initgap;
       
end

%Summarise results
fprintf("Final estimate: xk = (%.4f, %.4f, %.4f, %.4f) \n", x(k,1), x(k,2), x(k,3), x(k,4));
fprintf("True minimiser: x* = ");
disp(true_min);
fprintf(" \n");
fprintf("Condition number of A: %.4f\n", r);
fprintf("Total iterations performed: %d\n", k-1);
fprintf("Convergence rate for %d iterations: %d\n", k-1, delta^(2*(k-1)));

%Value of xk for each iteration
fprintf("xk for each iteration \n");
fprintf("k \n");
for i=1:1:k

    fprintf("%d| ", i-1);
    disp(x(i,:));

end

fprintf("-------------------------------\n");

%Convergence bounds of each iteration
fprintf("Convergence bounds for each iteration \n");
fprintf("k        convergence_rate     bounds\n");

for i=1:1:k

    fprintf("%d|           ", i-1);
    %fprintf("%d            %d\n", convk(i,1), convk(i,2));
    disp(convk(i,:));

end

disp(initgap);

end

