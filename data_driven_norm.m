function f = data_driven_norm(x)
    uc = 0.1:0.1:5.1;
    ur = 5.1:0.1:100.0;
    Hu = hankel(uc, ur);
    H = [Hu];
    w = [[15.1:0.1:20.1]'];
    f = norm(H*x'-w);
end
