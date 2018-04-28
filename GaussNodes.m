function [Xk,Ak] = GaussNodes(n)
%% Initialization
    n = n + 1;
    L = zeros(n+1,n+1);
    y = zeros(1,n);
    Xk = zeros(n-1,n-1);
    Ak = zeros(n-1,n-1);
    N = 1;
    syms x y;

%% Laguerre polynomial iteration
    L(1,n+1) = 1;
    L(2,n:n+1) = [-1,1];
    for i = 2:n
        d(1,n:n+1) = [-1, 1 + 2*(i-1)];
        dL = conv(d,L(i,:));
        dL = dL(end-n:end);
        L(i+1,:) = dL - (i-1)^2 * L(i-1,:);
    end

%% Solution
    for i = 2:n+1
        y(i-1) = 0;
        for j = 0:n
            y(i-1) = y(i-1) + L(i,n+1-j)*x^j;
        end
        Xk(i-1,1:i-1) = double(solve(y(i-1),x));
    end
    Xk = real(Xk);

%% Integral coefficient
    for i = 1:n-1
        N = N*i;
        for j = 1:i
            Lk = double(subs(y(i+1),x,Xk(i,j)));
            Ak(i,j) = N^2/Lk^2*Xk(i,j);
        end
    end
Ak = Ak(n-1,:);
Xk = Xk(n-1,1:n-1);
end

