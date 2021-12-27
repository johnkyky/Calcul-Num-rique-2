function [y] = spmv(AA, JA, IA, x)
    n = size(IA, 1)
    y = zeros(n - 1, 1)
    for i = 1 : n - 1
        for j = IA(i) : IA(i + 1) - 1
            y(i) = y(i) + AA(j) * x(JA(j))
        end
    end
end


function [A] = generate_sparse_random_matrix(m, n, coef)
    A = full(sprand(m, n, coef))
end


function [AA, JA, IA] = geToGb(A)
    B = sparse(A)
    [IJ, AA, nm]=spget(B)
    JA = IJ(:, 2)
    IA_tmp = []
    count = 0
    index = 1
    i = 1
    while i <= size(IJ, 1)
        if IJ(i, 1) == index
            count = count + 1
            i = i + 1
       else
            IA_tmp($ + 1) = count
            index = index + 1
            count = 0
        end
    end
    IA_tmp($ + 1) = count

    IA = []
    IA($ + 1) = 1
    for i = 2 : size(IA_tmp, 1)
        IA(i) = IA(i - 1) + IA_tmp(i - 1)
    end
    IA($ + 1) = IA($) + IA_tmp($)
end
