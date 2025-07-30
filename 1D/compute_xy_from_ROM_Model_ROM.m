function result = compute_xy_from_ROM_Model_ROM(Model, Opts, Xi)

N_xi = size(Xi,1);

result.xi = Xi;
result.U = sparse(N_xi, Opts.para.rus*Opts.para.n_f);
result.P = sparse(N_xi, Opts.para.rp*Opts.para.n_f);
result.Yu = sparse(N_xi, Opts.para.n_ou*Opts.para.n_f);
result.Yp = sparse(N_xi, Opts.para.n_op*Opts.para.n_f);

for i = 1:N_xi
    Opts.para.int_para.xi = Xi(i);
    [A_hat, B_hat, C_hat] = compute_ABC_from_Model_ROM(Model, Opts.para, Opts.ROM_coe_fun);
    x_hat = full(A_hat\B_hat);
    u = x_hat(1:Opts.para.rus,:);
    p = x_hat(Opts.para.rus+1:end,:);
    yu = C_hat{1} * u;
    yp = C_hat{2} * p;

    result.U(i,:) = u(:)';
    result.P(i,:) = p(:)';

    result.Yu(i,:) = yu(:)';
    result.Yp(i,:) = yp(:)';
end
end
