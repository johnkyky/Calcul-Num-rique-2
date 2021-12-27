getd('.')


x = []
y_ge = []
y_gb = []



SAMPLES = 10
index = 1
i = 10
while i <= 300
    disp(i)
    // initialisation variable
    A = generate_sparse_random_matrix(i, i, 0.5)
    [AA, JA, IA] = geToGb(A)
    v = rand(i, 1)
    
    // mesure de temps
    tic()
    for j = 1 : SAMPLES
        y = spmv(AA, JA, IA, v)
    end
    y_gb(index) = toc() / SAMPLES
    
    tic()
    for j = 1 : SAMPLES
        y_ = A * v
    end
    y_ge(index) = toc() / SAMPLES

    // validation mesure
    if abs(norm(y) - norm(y_)) > 1.0D-12
        disp("La ultiplication n est pas valide")
    end
    
    x(index) = i
    i = i + 10
    index = index + 1
end


clf()
xlabel("Taille matrice (n * n)")
ylabel("temps de calcul (seconde)")
plot("nl", x, [y_ge, y_gb])
legend(["stockage général", "stockage CSR"])


//Utilisation mémoire général vs CSR

/*
for i = 3 : 500
    x(i) = i
    y_ge(i) = i * i
    y_gb(i) = 2 * (i * 0.5) + i + 1
end

clf()
xlabel("Taille matrice (n * n)")
ylabel("nombre éléments")
plot("nl", x, [y_ge, y_gb])
legend(["stockage général", "stockage CSR"])

disp(y_ge($))
disp(y_gb($))
*/
