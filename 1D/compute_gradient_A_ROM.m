function gradient = compute_gradient_A_ROM(Opts)

int_coe = 2*Opts.Int.coe;
N_int_point = size(Opts.Int.point,1);

temp_A = cell(Opts.para.num.A, 1);
for j = 1:Opts.para.num.A
    temp_A{j} = zeros(Opts.para.r, Opts.para.r);
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

    for k = 1:Opts.para.num.A
        coe_fun_val = feval(Opts.ROM_coe_fun.A{k}, Opts.para.int_para);
        temp_A{k} = temp_A{k} - coe_fun_val * (x_hat_dual * temp_y * x_hat');
    end
end


gradient = [];
for j = 1:Opts.para.num.A
    temp_A{j} = (temp_A{j} + temp_A{j}')/2;
    temp_A{j} = temp_A{j}(tril(true(size(temp_A{j}))));
    gradient = [gradient; temp_A{j}(1:end-Opts.para.del)];
end

gradient = sparse(int_coe * gradient);

end