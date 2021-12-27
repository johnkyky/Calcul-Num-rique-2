n = 3
lim = 10e-6

A = toeplitz([2, -1, 0])
b = rand(n, 1)

eigens = spec(A)
disp("alpha optimum = ", 2 / (min(eigens) + max(eigens)))

residu_x1 = []
residu_x2 = []
residu_x3 = []
residu_x4 = []
residu_x5 = []
residu_x6 = []

residu_y1 = []
residu_y2 = []
residu_y3 = []
residu_y4 = []
residu_y5 = []
residu_y6 = []

x = zeros(n, 1)
i = 1
is_not_finish = 1
while is_not_finish
    res = b - A * x
    residu_y1(i) = norm(res)
    residu_x1(i) = i

    if(residu_y1(i) <= lim || i > 200)
        is_not_finish = 0
    end

    x = x + 0.1 * res
    i = i + 1
end


x = zeros(n, 1)
i = 1
is_not_finish = 1
while is_not_finish
    res = b - A * x
    residu_y2(i) = norm(res)
    residu_x2(i) = i

    if(residu_y2(i) <= lim || i > 200)
        is_not_finish = 0
    end

    x = x + 0.2 * res
    i = i + 1
end


x = zeros(n, 1)
i = 1
is_not_finish = 1
while is_not_finish
    res = b - A * x
    residu_y3(i) = norm(res)
    residu_x3(i) = i

    if(residu_y3(i) <= lim || i > 200)
        is_not_finish = 0
    end

    x = x + 0.3 * res
    i = i + 1
end

x = zeros(n, 1)
i = 1
is_not_finish = 1
while is_not_finish
    res = b - A * x
    residu_y4(i) = norm(res)
    residu_x4(i) = i

    if(residu_y4(i) <= lim || i > 200)
        is_not_finish = 0
    end

    x = x + 0.4 * res
    i = i + 1
end


x = zeros(n, 1)
i = 1
is_not_finish = 1
while is_not_finish
    res = b - A * x
    residu_y5(i) = norm(res)
    residu_x5(i) = i

    if(residu_y5(i) <= lim || i > 200)
        is_not_finish = 0
    end

    x = x + 0.5 * res
    i = i + 1
end



x = zeros(n, 1)
i = 1
is_not_finish = 1
while is_not_finish
    res = b - A * x
    residu_y6(i) = norm(res)
    residu_x6(i) = i

    if(residu_y6(i) <= lim || i > 200)
        is_not_finish = 0
    end

    x = x + 0.6 * res
    i = i + 1
end



//Residu limit
size_max = max(size(residu_y1, 1), size(residu_y2, 1), size(residu_y3, 1), size(residu_y4, 1), size(residu_y5, 1), size(residu_y6, 1))
residu_xlim = [1, size_max]
residu_ylim = [lim, lim]


clf()
xlabel("nombre d itération")
ylabel("norme du residu")
plot("nl", residu_x1, residu_y1, residu_x2, residu_y2, residu_x3, residu_y3, residu_x4, residu_y4, residu_x5, residu_y5, residu_x6, residu_y6, residu_xlim, residu_ylim)
legend(["alpha = 0.1", "alpha = 0.2", "alpha = 0.3", "alpha = 0.4", "alpha = 0.5", "alpha = 0.6", "erreur limite : 10e-6"])


