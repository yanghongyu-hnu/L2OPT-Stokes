function gradient = compute_gradient_B_ROM(Opts)

int_coe = 2 * Opts.Int.coe;
N_int_point = size(Opts.Int.point,1);

temp_B = cell(Opts.para.num.B, 1);
for j = 1:Opts.para.num.B
    temp_B{j} = zeros(Opts.para.r, Opts.para.n_f);
end

for i = 1:N_int_point
    Opts.para.int_para.xi = Opts.Int.point(i);
    [A_hat, B_hat, C_hat] = compute_ABC_from_Model_ROM(Opts.iter.ROM_Model, Opts.para, Opts.ROM_coe_fun);
    blk_C_hat = blkdiag(C_hat{1},C_hat{2});

    x_hat_dual = full(A_hat'\blk_C_hat');
    x_hat = full(A_hat\B_hat);
    y_hat = full(blk_C_hat * x_hat);

    y_current_int_point = full([reshape(Opts.Int.Y_int_point.Yu(i,:), Opts.para.n_ou, Opts.para.n_f); reshape(Opts.Int.Y_int_point.Yp(i,:), Opts.para.n_op, Opts.para.n_f)]);

    temp_y = y_hat - y_current_int_point;

    for k = 1:Opts.para.num.B
        coe_fun_val = feval(Opts.ROM_coe_fun.B{k}, Opts.para.int_para);
        temp_B{k} = temp_B{k} + coe_fun_val * (x_hat_dual * temp_y);
    end

end

gradient = [];
for j = 1:Opts.para.num.B
    gradient = [gradient; temp_B{j}(:)];
end

gradient = sparse(int_coe * gradient);

end