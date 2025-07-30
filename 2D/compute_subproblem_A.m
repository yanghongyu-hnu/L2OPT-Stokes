function [ofun_val_1, current_Model, x1, step_size, gradient_diff_A, x_diff_A, L_BFGS_A, flag_1, i_inner] = compute_subproblem_A(max_rate, ofun_val_0, L_BFGS_A,  gradient_0_diff_A, x_0_diff_A, iter, Opts)
max_iter = 200;

x0 = Opts.iter.x;
current_Model = Opts.iter.ROM_Model;
x0_mat = mat2cell(x0, Opts.para.partition, 1);
x_A_0 = x0_mat{1};
inexact = max(1,sqrt(s_inner_product(x_A_0, x_A_0, Opts))) * 1e-13;
gradient_A_0 = compute_gradient_A_ROM(Opts);

if max_rate >= 1e-4 || iter <= Opts.L2_Opts.m || L_BFGS_A.iter_m <= Opts.L2_Opts.m || isempty(x_0_diff_A)
    d_A = -gradient_A_0;
    flag_1 = 0;
else
    gamma = (s_inner_product(x_0_diff_A, gradient_0_diff_A, Opts)/s_inner_product(gradient_0_diff_A, gradient_0_diff_A, Opts));
    d_A = -L_bfgs_A(L_BFGS_A.rho, gradient_A_0, L_BFGS_A.S, L_BFGS_A.Y, gamma, Opts);
    flag_1 = 1;
end

x_init = x_A_0;
x_temp_0 = [x_A_0;x0_mat{2};x0_mat{3}];

for i_inner = 1:max_iter
    d_A = d_A/sqrt(s_inner_product(d_A, d_A, Opts));
    alpha = Opts.L2_Opts.init_step_size;
    x_A_1 = x_A_0 + alpha * d_A;
    x_temp_1 = [x_A_1;x0_mat{2};x0_mat{3}];
    temp_gradient = Opts.L2_Opts.Armijo.sigma * s_inner_product(gradient_A_0, d_A, Opts);
    Opts.iter.ROM_Model = convert_vector_to_Model_ROM(Opts.para, x_temp_1);
    ofun_val_1 = compute_ofun_ROM(Opts);

    while ofun_val_1.y > ofun_val_0.y + alpha * temp_gradient
        alpha = alpha * Opts.L2_Opts.Armijo.rho;
        x_A_1 = x_A_0 + alpha * d_A;
        x_temp_1 = [x_A_1;x0_mat{2};x0_mat{3}];
        Opts.iter.ROM_Model = convert_vector_to_Model_ROM(Opts.para, x_temp_1);
        ofun_val_1 = compute_ofun_ROM(Opts);
        if alpha <= 1e-14
            break;
        end
    end

    if ofun_val_1.y < ofun_val_0.y
        flag_2 = 0;
        x_diff_A = x_A_1 - x_A_0;
        gradient_A_1 = compute_gradient_A_ROM(Opts) + x_A_1 - x_init;
        gradient_diff_A = gradient_A_1 - gradient_A_0;
        temp = s_inner_product(gradient_diff_A, x_diff_A, Opts);

        if temp > 0
            if L_BFGS_A.iter_m <= Opts.L2_Opts.m
                L_BFGS_A.Y(:,L_BFGS_A.iter_m) = gradient_diff_A;
                L_BFGS_A.S(:,L_BFGS_A.iter_m) = x_diff_A;
                L_BFGS_A.rho(L_BFGS_A.iter_m) = 1/temp;
                L_BFGS_A.iter_m = L_BFGS_A.iter_m + 1;
            else
                L_BFGS_A.rho(1:Opts.L2_Opts.m-1) = L_BFGS_A.rho(2:Opts.L2_Opts.m); L_BFGS_A.rho(Opts.L2_Opts.m) = 1/temp;
                L_BFGS_A.S(:,1:Opts.L2_Opts.m-1) = L_BFGS_A.S(:,2:Opts.L2_Opts.m); L_BFGS_A.S(:,Opts.L2_Opts.m) = x_diff_A;
                L_BFGS_A.Y(:,1:Opts.L2_Opts.m-1) = L_BFGS_A.Y(:,2:Opts.L2_Opts.m); L_BFGS_A.Y(:,Opts.L2_Opts.m) = gradient_diff_A;
            end
        end

        if sqrt(s_inner_product(x_diff_A, x_diff_A, Opts)) < inexact || i_inner == max_iter
            x1 = x_temp_1;
            current_Model = Opts.iter.ROM_Model;
            flag_1 = 1;
            break;
        end

        if L_BFGS_A.iter_m > Opts.L2_Opts.m
            flag_1 = 1;
            gamma = (s_inner_product(x_diff_A, gradient_diff_A, Opts)/s_inner_product(gradient_diff_A, gradient_diff_A, Opts));
            d_A = -L_bfgs_A(L_BFGS_A.rho, gradient_A_1, L_BFGS_A.S, L_BFGS_A.Y, gamma, Opts);
        else
            flag_1 = 0;
            d_A = -gradient_A_1;
        end
        gradient_A_0 = gradient_A_1;
        x_A_0 = x_A_1;
        ofun_val_0 = ofun_val_1;
        x_temp_0 = x_temp_1;
    else
        if i_inner == 1
            if flag_1 == 0 && L_BFGS_A.iter_m > Opts.L2_Opts.m
                if isempty(x_0_diff_A)
                    flag_1 = 0;
                    x1 = x0;
                    ofun_val_1 = ofun_val_0;
                    gradient_diff_A = [];
                    x_diff_A = [];
                    break
                end
                flag_1 = 1;
                gamma = (s_inner_product(x_0_diff_A, gradient_0_diff_A, Opts)/s_inner_product(gradient_0_diff_A, gradient_0_diff_A, Opts));
                d_A = -L_bfgs_A(L_BFGS_A.rho, gradient_A_0, L_BFGS_A.S, L_BFGS_A.Y, gamma, Opts);
            else
                flag_1 = 0;
                d_A = -gradient_A_0;
            end
            flag_2 = 1;
        else
            if flag_2 == 0
                if flag_1 == 0 && L_BFGS_A.iter_m > Opts.L2_Opts.m
                    flag_1 = 1;
                    gamma = (s_inner_product(x_diff_A, gradient_diff_A, Opts)/s_inner_product(gradient_diff_A, gradient_diff_A, Opts));
                    d_A = -L_bfgs_A(L_BFGS_A.rho, gradient_A_0, L_BFGS_A.S, L_BFGS_A.Y, gamma, Opts);
                else
                    flag_1 = 0;
                    d_A = -gradient_A_0;
                end
                flag_2 = 1;
            else
                flag_1 = 0;
                x1 = x_temp_0;
                ofun_val_1 = ofun_val_0;
                current_Model = convert_vector_to_Model_ROM(Opts.para, x1);
                if i_inner <= 2
                    gradient_diff_A = [];
                    x_diff_A = [];
                end
                break;
            end
        end
    end
end
if i_inner == max_iter
    x1 = x_temp_1;
    current_Model = convert_vector_to_Model_ROM(Opts.para, x1);
    flag_1 = 2;
end
step_size = alpha;
end