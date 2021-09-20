%%
% #INPUT MATRIX

matA = magic(9);
%%
% #DETERMINANT

% METHOD_1:Cofactor
D_1 = 0;

for n = 1:size(matA, 2) %Iteration for nth column in 1st row
    D_1 = D_1 + cof(matA, n) * matA(1, n); % D=sum(cofactor*entries)
end

D_1
% METHOD_2:Product of diagonal entries in row echelon form(REF)
E = GE(matA)
D_2 = prod(diag(E))
%%
% CHECK

det(matA)
round(det(matA) - D_1, -2) == 0 % METHOD_1: 100의 자리까지 일치
det(matA) - D_2 == 0 % METHOD_2: 모두 일치
%%
% #INVERSE

% METHOD: Gauss-Jordan
augA = [matA eye(9)]
RE = GJE(augA)
invA = RE(:, 10:end)
%%
% CHECK

inv(matA)
round(invA - inv(matA), 14) == 0 % GJE: 1.0e-14의 자리까지 일치
%%
% FUNCTIONS

function M = minor(A, col) % Minor of a=A(1,col)
    M = 0;
    A(1, :) = [];
    A(:, col) = [];

    if size(A) ~= 2

        for n = 1:size(A, 2)
            a = A(1, n);
            M = M + (-1)^(1 + n) * a * minor(A, n);
        end

    else
        M = A(1) * A(4) - A(2) * A(3); % 2x2 determinant
    end

end

function C = cof(A, col) % Cofactor of a=A(1,col)
    C = (-1)^(1 + col) * minor(A, col);
end

function E = GE(A) % Gaussian Elimination

    for n = 1:size(A, 2)

        for m = n + 1:size(A, 1)
            A(m, :) = A(m, :) - A(n, :) * A(m, n) / A(n, n);
            E = A;
        end

    end

end

function RE = GJE(A) % Gauss-Jordan Elimination

    for n = 1:size(A, 1)

        for m = [n + 1:size(A, 1) 1:n - 1]
            A(m, :) = A(m, :) - A(n, :) * A(m, n) / A(n, n);
        end

        A(n, :) = A(n, :) / A(n, n); % Scaling
        RE = A;
    end

end
