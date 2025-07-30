function Opts = default_L2_Opt_ROM(Opts)

if Opts.para.num.A == size(Opts.initial_point.A, 1)
    initial_point.A = Opts.initial_point.A;
end

if Opts.para.num.B == size(Opts.initial_point.B, 1)
    initial_point.B = Opts.initial_point.B;
end

if Opts.para.num.C == size(Opts.initial_point.C, 1)
    initial_point.C = Opts.initial_point.C;
end

initial_point.blk = Opts.blk;
Opts.iter.ROM_Model = initial_point;

Opts.int_init_flag = 1;
Opts.Int.N_int_point = size(Opts.Int.point,1)^(size(Opts.Int.point,2));
Opts.Int.coe = 1/Opts.Int.N_int_point;

Opts.para.r = Opts.para.rus + Opts.para.rp;
Opts.para.del = (Opts.para.rp*(Opts.para.rp+1))/2;

Opts.para.size_A = [Opts.para.r; Opts.para.r*(Opts.para.r+1)/2 - Opts.para.del];
Opts.para.size_B = Opts.para.r * Opts.para.n_f;
Opts.para.size_C = Opts.para.n_ou * Opts.para.rus + Opts.para.n_op * Opts.para.rp;

Opts.para.partition = [Opts.para.size_A(2) * Opts.para.num.A, Opts.para.num.B * Opts.para.size_B, Opts.para.num.C * Opts.para.size_C];

Opts.para.total_num_var = sum(Opts.para.partition);

Opts.Int.Y_int_point.Yu = Opts.train_data.Yu_modified;
Opts.Int.Y_int_point.Yp = Opts.train_data.Yp;

n = Opts.para.size_A(1,1);
d = ones(n, 1);
temp = 2 * ones(n) - spdiags(d, 0, n, n);
temp = temp(tril(true(size(temp))));
temp = temp(1:end - Opts.para.del);
temp = temp(:, ones(1, Opts.para.num.A));
Opts.para_vector = temp(:);


end