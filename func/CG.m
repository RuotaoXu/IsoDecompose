function [ X ] = CG(A, B, X0, n_loop, tol, P, valid_func)
%% Conjugate Gradient method
% Version: 1.3
%% Parameter Setting
if ~exist('X0','var') || isempty(X0)
    X0 = zeros(size(A,2),size(B,1));
end
if ~exist('P','var')
    P = [];
end
if ~exist('n_loop','var') || isempty(n_loop)
    n_loop = 100;
end
if ~exist('tol','var') || isempty(tol)
    tol = realmin;
end

%% Iteration
if isempty(P)
    r0 = A(X0) - B;
    Q = -r0;
    X = X0;
    
    if norm(r0(:),2)/numel(r0)<tol return; end
    for iter = 1:n_loop
        AQ = A(Q);
        r00 = r0(:)'*r0(:);
        alpha = r00/(Q(:)'*AQ(:));
        X = X + alpha * Q;
        r1 = r0 + alpha * AQ;
        r11 = r1(:)'*r1(:);
%         r=A(X) - B;
%         fprintf('%.2f\n',r(:)'*r(:));
        if r11< tol || (exist('valid_func','var') && ~isempty(valid_func) && valid_func(X))
            break;
        end
        beta = r11 /r00;
        Q = -r1 + beta * Q;
        r0 = r1;
    end
else
    R = B - A(X0);
    Z = P(R);
    Q = Z;
    X = X0;
    for iter = 1: n_loop
        AQ = A(Q);
        zr = sums(R.*Z);
        alpha = zr/sums(Q.*AQ);
        X = X + alpha*Q;
        R = R - alpha*AQ;
        if norm(R,'fro') < tol || (exist('valid_func','var') && ~isempty(valid_func) && valid_func(X))
            break;
        end
        Z = P(R);
        beta = sums(Z.*R)/zr;
        Q = Z + beta * Q;
    end
end

end