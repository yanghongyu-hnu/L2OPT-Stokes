function ofun_val = compute_ofun_ROM(Opts)

temp_yu = 0;
temp_yp = 0;

N_int_point = size(Opts.Int.point,1);

for i = 1:N_int_point
    start = (i-1) * N_int_point;
    for j = 1:N_int_point
        row = start + j;
        Opts.para.int_para.xi = [Opts.Int.point(i,1); Opts.Int.point(j,2)];
        [A_hat, B_hat, C_hat] = compute_ABC_from_Model_ROM(Opts.iter.ROM_Model, Opts.para, Opts.ROM_coe_fun);

        x_hat = full(A_hat\B_hat);

        u = x_hat(1:Opts.para.rus,:);
        p = x_hat(Opts.para.rus+1:end,:);

        yu_hat = C_hat{1} * u;
        yp_hat = C_hat{2} * p;

        yu = full(reshape(Opts.Int.Y_int_point.Yu(row,:), Opts.para.n_ou, Opts.para.n_f));
        yp = full(reshape(Opts.Int.Y_int_point.Yp(row,:), Opts.para.n_op, Opts.para.n_f));

        temp_yu = temp_yu + sum(sum((yu - yu_hat).^2));
        temp_yp = temp_yp + sum(sum((yp - yp_hat).^2));
    end
end

ofun_val.yp = Opts.Int.coe * temp_yp;
ofun_val.yu = Opts.Int.coe * temp_yu;
ofun_val.y = ofun_val.yp + ofun_val.yu;

end