getd('.')


x = []
y_ldlt = []
y_lu = []


SAMPLES = 10
index = 1
i = 10
while i <= 100
    disp(i)
    // initialisation variable
    A = rand(i, i)
    A = symetrise(A)
    
    // mesure de temps
    tic()
    for j = 1 : SAMPLES
        [L, D] = myldlt(A)
    end
    y_ldlt(index) = toc() / SAMPLES

    // validation mesure
    if abs(norm(A) - norm(L * D * L')) > 1.0D-14
        disp("La factorisation LDLT n est pas valide")
    end

     // mesure de temps
    tic()
    for j = 1 : SAMPLES
        [L, U] = mylu_opti(A)
    end
    y_lu(index) = toc() / SAMPLES


    x(index) = i
    index = index + 1
    i = i + 10
end

clf()
xlabel("Taille matrice (n * n)")
ylabel("temps de calcul (seconde)")
plot(x, [y_ldlt, y_lu])
legend(["LDLt", "LU"])
