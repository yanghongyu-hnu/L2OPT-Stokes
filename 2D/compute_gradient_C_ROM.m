function gradient = compute_gradient_C_ROM(Opts)

int_coe = 2 * Opts.Int.coe;
N_int_point = size(Opts.Int.point,1);

temp_C = cell(Opts.para.num.C, Opts.para.num_output);
for j = 1:Opts.para.num.C
    temp_C{j,1} = zeros(Opts.para.n_ou, Opts.para.rus);
    temp_C{j,2} = zeros(Opts.para.n_op, Opts.para.rp);
end

for i = 1:N_int_point
    start = (i-1) * N_int_point;
    for j = 1:N_int_point
        row = start + j;
        Opts.para.int_para.xi = [Opts.Int.point(i,1); Opts.Int.point(j,2)];
        [A_hat, B_hat, C_hat] = compute_ABC_from_Model_ROM(Opts.iter.ROM_Model, Opts.para, Opts.ROM_coe_fun);
        blk_C_hat = blkdiag(C_hat{1},C_hat{2});

        x_hat = full(A_hat\B_hat);
        y_hat = full(blk_C_hat * x_hat);

        y_current_int_point = full([reshape(Opts.Int.Y_int_point.Yu(row,:), Opts.para.n_ou, Opts.para.n_f); reshape(Opts.Int.Y_int_point.Yp(row,:), Opts.para.n_op, Opts.para.n_f)]);

        temp_y = y_hat - y_current_int_point;
        temp_yu = temp_y(1:Opts.para.n_ou,:);
        temp_yp = temp_y(Opts.para.n_ou+1:end,:);

        for k = 1:Opts.para.num.C
            coe_fun_val = feval(Opts.ROM_coe_fun.C{k}, Opts.para.int_para);
            temp_C{k,1} = temp_C{k,1} + coe_fun_val * (temp_yu * x_hat(1:Opts.para.rus)');
            temp_C{k,2} = temp_C{k,2} + coe_fun_val * (temp_yp * x_hat(Opts.para.rus+1:end)');
        end
    end
end

gradient = [];
for j = 1:Opts.para.num.C
    gradient = [gradient; temp_C{j,1}(:)];
    gradient = [gradient; temp_C{j,2}(:)];
end

gradient = sparse(int_coe * gradient);

end