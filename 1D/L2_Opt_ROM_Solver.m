function result = L2_Opt_ROM_Solver(Opts)

tic;
tol = clock;

Opts = default_L2_Opt_ROM(Opts);

y_hat_0 = compute_y_from_current_Model_ROM(Opts.iter.ROM_Model, Opts.para, Opts.ROM_coe_fun, Opts.Int.point);
ofun_val_0 = compute_ofun_ROM(Opts);

fprintf('\n--------------------------------------------------------------');
fprintf('\n  iter | iter.A| iter.B| iter.C|  ||yu0-y1||  |   ||yu1||   |    rate    |   ||yp0-yp1||   |   ||yp1||   |    rate    |   ||y0-y1||   |   ||y1||   |    rate    |             ofun_val             | step_size  | step_size  | step_size  | time  | total time  |\n');
fprintf(' %3d   |   -   |   -   |   -   |       -      |      -      |     -      |         -       |      -      |     -      |        -      |      -     |     -      | %2.4e|%2.4e|%2.4e |     -      |     -      |     -      |   -   |    %s    |\n', 0, ofun_val_0.y, ofun_val_0.yu, ofun_val_0.yp, num2str(etime(clock,tol)));

x_0 = convert_Model_to_vector_ROM(Opts.iter.ROM_Model, Opts.para);
Opts.iter.x = x_0;

L_BFGS.A.rho = zeros(1,Opts.L2_Opts.m);
L_BFGS.A.S = zeros(Opts.para.partition(1), Opts.L2_Opts.m);
L_BFGS.A.Y = zeros(Opts.para.partition(1), Opts.L2_Opts.m);
L_BFGS.A.iter_m = 1;

L_BFGS.B.rho = zeros(1,Opts.L2_Opts.m);
L_BFGS.B.S = zeros(Opts.para.partition(2), Opts.L2_Opts.m);
L_BFGS.B.Y = zeros(Opts.para.partition(2), Opts.L2_Opts.m);
L_BFGS.B.iter_m = 1;

L_BFGS.C.rho = zeros(1,Opts.L2_Opts.m);
L_BFGS.C.S = zeros(Opts.para.partition(3), Opts.L2_Opts.m);
L_BFGS.C.Y = zeros(Opts.para.partition(3), Opts.L2_Opts.m);
L_BFGS.C.iter_m = 1;

gradient_diff.A = [];
gradient_diff.B = [];
gradient_diff.C = [];
x_diff.A = [];
x_diff.B = [];
x_diff.C = [];

max_rate = 1;

for i = 1:Opts.L2_Opts.maxit
    t_i = clock;

    [ofun_val_1_A, Opts.iter.ROM_Model, Opts.iter.x, step_size.A, gradient_diff.A, x_diff.A, L_BFGS.A, flag.A, inner_iter.A] = compute_subproblem_A(max_rate, ofun_val_0, L_BFGS.A, gradient_diff.A, x_diff.A, i, Opts);
    [ofun_val_1_B, Opts.iter.ROM_Model, Opts.iter.x, step_size.B, gradient_diff.B, x_diff.B, L_BFGS.B, flag.B, inner_iter.B] = compute_subproblem_B(max_rate, ofun_val_1_A, L_BFGS.B, gradient_diff.B, x_diff.B, i, Opts);
    [ofun_val_1, Opts.iter.ROM_Model, x_1, step_size.C, gradient_diff.C, x_diff.C, L_BFGS.C, flag.C, inner_iter.C] = compute_subproblem_C(max_rate, ofun_val_1_B, L_BFGS.C, gradient_diff.C, x_diff.C, i, Opts);

    y_hat_1 = compute_y_from_current_Model_ROM(Opts.iter.ROM_Model, Opts.para, Opts.ROM_coe_fun, Opts.Int.point);

    temp_1.norm_yu = sqrt(sum(sum(full(((y_hat_0.Yu_hat - y_hat_1.Yu_hat)').^2))) * Opts.Int.coe);
    temp_2.norm_yu = sqrt(sum(sum(full((y_hat_1.Yu_hat').^2))) * Opts.Int.coe);
    rate.norm_yu = temp_1.norm_yu/(1 + temp_2.norm_yu);

    temp_1.norm_yp = sqrt(sum(sum(full(((y_hat_0.Yp_hat - y_hat_1.Yp_hat)').^2))) * Opts.Int.coe);
    temp_2.norm_yp = sqrt(sum(sum(full((y_hat_1.Yp_hat').^2))) * Opts.Int.coe);
    rate.norm_yp = temp_1.norm_yp/(1 + temp_2.norm_yp);

    temp_1.norm_y = sqrt(temp_1.norm_yu^2 + temp_1.norm_yp^2);
    temp_2.norm_y = sqrt(temp_2.norm_yu^2 + temp_2.norm_yp^2);
    rate.norm_y = temp_1.norm_y/(1 + temp_2.norm_y);

    fprintf(' %3d   | %3d   | %3d   | %3d   |  %2.4e  | %2.4e  | %2.4e |    %2.4e   | %2.4e  | %2.4e |   %2.4e  | %2.4e | %2.4e | %2.4e|%2.4e|%2.4e | %2.4e | %2.4e | %2.4e | %s |    %s    |\n', i, inner_iter.A, inner_iter.B, inner_iter.C, temp_1.norm_yu, temp_2.norm_yu, rate.norm_yu, temp_1.norm_yp, temp_2.norm_yp, rate.norm_yp, temp_1.norm_y, temp_2.norm_y, rate.norm_y, ofun_val_1.y, ofun_val_1.yu, ofun_val_1.yp, step_size.A, step_size.B, step_size.C, num2str(etime(clock,t_i)), num2str(etime(clock,tol)));

    max_rate = max([rate.norm_y,rate.norm_yp,rate.norm_yu]);
    min_rate = min([rate.norm_y,rate.norm_yp,rate.norm_yu]);
    max_val = max([temp_1.norm_y,temp_1.norm_yu,temp_1.norm_yp]);
    if max_rate <= Opts.L2_Opts.tol
        fprintf('--------------------------------------------------------------\n');
        fprintf('rate_norm_y = %2.4e\n',rate.norm_y);
        fprintf('rate_norm_yu = %2.4e\n',rate.norm_yu);
        fprintf('rate_norm_yp = %2.4e\n',rate.norm_yp);
        fprintf('max_rate = %2.4e < %2.4e \n',max_rate,Opts.L2_Opts.tol);
        fprintf('min_rate = %2.4e\n',min_rate);
        fprintf('max_val = %2.4e\n', max_val);
        fprintf('用时：%ss\n',num2str(etime(clock,tol)));
        fprintf('--------------------------------------------------------------\n');
        break;
    end

    Opts.iter.x = x_1;
    Opts.iter.ROM_Model = convert_vector_to_Model_ROM(Opts.para, x_1);
    y_hat_0 = y_hat_1;
    ofun_val_0 = ofun_val_1;
end

result = recover_from(Opts.iter.ROM_Model, Opts.para);


end