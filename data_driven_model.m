function data_driven_model
    A = zeros(950);
    b = ones(950,1);
    [a ,f] = fmincon(@data_driven_norm, zeros(1, 950), A, b);
    xc = 2*(0.1:0.1:5.1)-3;
    xr = 2*(5.1:0.1:100)-3;
    Hx = hankel(xc, xr);
    y = Hx * a';
    plot(y);
end