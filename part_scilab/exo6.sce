residu_xj = []
residu_j = []
residu_xgs = []
residu_gs = []
x_t = []
y_t_j = []
y_t_gs = []

lim = 10e-12


//////////////////////////////////////////
//       Mesure nombre iteration        //
//////////////////////////////////////////


n = 3

vec = zeros(1, n)
vec(1) = 2
vec(2) = -1
A = toeplitz(vec)
b = rand(n, 1)


//jacobi
Mj_inv = inv(diag(diag(A)))
xj = zeros(n, 1)
//Gauss Seidel
Mgs_inv = inv(tril(A))
xgs = zeros(n, 1)

//jacobi
i = 1
is_not_finish = 1
while is_not_finish
    res_j = b - A * xj
    residu_j(i) = norm(res_j)
    residu_xj(i) = i

    if(residu_j(i) <= lim)
        is_not_finish = 0
    end

    xj = xj + Mj_inv * res_j
    i = i + 1
end


//Gauss Seidel
i = 1
is_not_finish = 1
tic()
while is_not_finish
    res_gs = b - A * xgs
    residu_gs(i) = norm(res_gs)
    residu_xgs(i) = i

    if(residu_gs(i) <= lim)
        is_not_finish = 0
    end

    xgs = xgs + Mgs_inv * res_gs
    i = i + 1
end
disp("temps gauss seidel", toc())


//Residu limit
size_max = max(size(residu_j, 1), size(residu_gs, 1))
residu_xlim = [1, size_max + 15]
residu_lim = [lim, lim]

clf()
xlabel("nombre d itération")
ylabel("norme du residu")
plot("nl", residu_xj, residu_j, residu_xgs, residu_gs, residu_xlim, residu_lim)
legend(["methode de Jacobi", "methode de Gauss Seidel", "erreur limit : 10e-12"])



//////////////////////////////////////////
//       Mesure temps d'exécution       //
//////////////////////////////////////////
/*
x_t = []
y_t_j = []
y_t_gs = []

n = 3
while n < 16
    vec = zeros(1, n)
    vec(1) = 2
    vec(2) = -1
    A = toeplitz(vec)
    b = rand(n, 1)

    //jacobi
    xj = zeros(n, 1)
    //Gauss Seidel
    Mgs_inv = inv(tril(A))
    xgs = zeros(n, 1)
    
    //jacobi
    i = 1
    is_not_finish = 1
    tic()
    for k = 1 : 10
        while is_not_finish
            res_j = b - A * xj
            residu_j(i) = norm(res_j)
    
            if(residu_j(i) <= lim)
                is_not_finish = 0
            end
    
            xj = xj + 0.5 * res_j
            i = i + 1
        end
    end
    y_t_j(n - 2) =  toc() / 10
    
    
    //gauss seidel
    i = 1
    is_not_finish = 1
    tic()
    for k = 1 : 10
        while is_not_finish
            res_gs = b - A * xj
            residu_gs(i) = norm(res_gs)
    
            if(residu_gs(i) <= lim)
                is_not_finish = 0
            end
    
            xgs = xgs + Mgs_inv * res_gs
            i = i + 1
        end
    end
    y_t_gs(n - 2) =  toc() / 10
    
    disp(n)
    x_t(n - 2) = n
    n = n + 1
end

clf()
xlabel("taille du problème")
ylabel("temps de résolution")
plot("nl", x_t, [y_t_j, y_t_gs])
legend(["methode de Jacobi", "methode de Gauss Seidel"])*/










