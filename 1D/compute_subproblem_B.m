function [ofun_val_1, current_Model, x1, step_size, gradient_diff_B, x_diff_B, L_BFGS_B, flag_1, i_inner] = compute_subproblem_B(max_rate, ofun_val_0, L_BFGS_B, gradient_0_diff_B, x_0_diff_B, iter, Opts)
max_iter = 100;

x0 = Opts.iter.x;
current_Model = Opts.iter.ROM_Model;
x0_mat = mat2cell(x0, Opts.para.partition, 1);
x_B_0 = x0_mat{2};
gradient_B_0 = compute_gradient_B_ROM(Opts);
inexact = 1/2 * max(1,norm(x_B_0)) * 1e-11;

if max_rate >= 1e-4 || L_BFGS_B.iter_m <= Opts.L2_Opts.m || isempty(x_0_diff_B)
    d_B = -gradient_B_0;
    flag_1 = 0;
else
    gamma = (x_0_diff_B'*gradient_0_diff_B)/(gradient_0_diff_B'*gradient_0_diff_B);
    d_B = -L_bfgs_B(L_BFGS_B.rho, gradient_B_0, L_BFGS_B.S, L_BFGS_B.Y, gamma, Opts);
    flag_1 = 1;
end

x_init = x_B_0;
x_temp_0 = [x0_mat{1};x_B_0;x0_mat{3}];

for i_inner = 1:max_iter
    d_B = d_B/norm(d_B);
    alpha = Opts.L2_Opts.init_step_size;
    x_B_1 = x_B_0 + alpha * d_B;
    x_temp_1 = [x0_mat{1};x_B_1;x0_mat{3}];
    temp_gradient = Opts.L2_Opts.Armijo.sigma * (gradient_B_0' * d_B);
    Opts.iter.ROM_Model = convert_vector_to_Model_ROM(Opts.para, x_temp_1);
    ofun_val_1 = compute_ofun_ROM(Opts);

    while ofun_val_1.y > ofun_val_0.y + alpha * temp_gradient
        alpha = alpha * Opts.L2_Opts.Armijo.rho;
        x_B_1 = x_B_0 + alpha * d_B;
        x_temp_1 = [x0_mat{1};x_B_1;x0_mat{3}];
        Opts.iter.ROM_Model = convert_vector_to_Model_ROM(Opts.para, x_temp_1);
        ofun_val_1 = compute_ofun_ROM(Opts);
        if alpha <= 1e-14
            break;
        end
    end

    if ofun_val_1.y < ofun_val_0.y
        flag_2 = 0;
        x_diff_B = x_B_1 - x_B_0;
        gradient_B_1 = compute_gradient_B_ROM(Opts) + x_B_1 - x_init;
        gradient_diff_B = gradient_B_1 - gradient_B_0;
        temp = gradient_diff_B'*x_diff_B;

        if temp > 0
            if L_BFGS_B.iter_m <= Opts.L2_Opts.m
                L_BFGS_B.Y(:,L_BFGS_B.iter_m) = gradient_diff_B;
                L_BFGS_B.S(:,L_BFGS_B.iter_m) = x_diff_B;
                L_BFGS_B.rho(L_BFGS_B.iter_m) = 1/temp;
                L_BFGS_B.iter_m = L_BFGS_B.iter_m + 1;
            else
                L_BFGS_B.rho(1:Opts.L2_Opts.m-1) = L_BFGS_B.rho(2:Opts.L2_Opts.m); L_BFGS_B.rho(Opts.L2_Opts.m) = 1/temp;
                L_BFGS_B.S(:,1:Opts.L2_Opts.m-1) = L_BFGS_B.S(:,2:Opts.L2_Opts.m); L_BFGS_B.S(:,Opts.L2_Opts.m) = x_diff_B;
                L_BFGS_B.Y(:,1:Opts.L2_Opts.m-1) = L_BFGS_B.Y(:,2:Opts.L2_Opts.m); L_BFGS_B.Y(:,Opts.L2_Opts.m) = gradient_diff_B;
            end
        end

        if norm(x_diff_B) < inexact || i_inner == max_iter
            x1 = x_temp_1;
            current_Model = Opts.iter.ROM_Model;
            flag_1 = 1;
            break;
        end

        if L_BFGS_B.iter_m > Opts.L2_Opts.m
            flag_1 = 1;
            gamma = (x_diff_B'*gradient_diff_B)/(gradient_diff_B'*gradient_diff_B);
            d_B = -L_bfgs_B(L_BFGS_B.rho, gradient_B_1, L_BFGS_B.S, L_BFGS_B.Y, gamma, Opts);
        else
            flag_1 = 0;
            d_B = -gradient_B_1;
        end
        gradient_B_0 = gradient_B_1;
        x_B_0 = x_B_1;
        ofun_val_0 = ofun_val_1;
        x_temp_0 = x_temp_1;
    else
        if i_inner == 1
            if flag_1 == 0 && L_BFGS_B.iter_m > Opts.L2_Opts.m
                if isempty(x_0_diff_B)
                    flag_1 = 0;
                    x1 = x0;
                    ofun_val_1 = ofun_val_0;
                    gradient_diff_B = [];
                    x_diff_B = [];
                    break;
                end
                flag_1 = 1;
                gamma = (x_0_diff_B'*gradient_0_diff_B)/(gradient_0_diff_B'*gradient_0_diff_B);
                d_B = -L_bfgs_B(L_BFGS_B.rho, gradient_B_0, L_BFGS_B.S, L_BFGS_B.Y, gamma, Opts);
            else
                flag_1 = 0;
                d_B = -gradient_B_0;
            end
            flag_2 = 1;
        else
            if flag_2 == 0
                if flag_1 == 0 && L_BFGS_B.iter_m > Opts.L2_Opts.m
                    flag_1 = 1;
                    gamma = (x_diff_B'*gradient_diff_B)/(gradient_diff_B'*gradient_diff_B);
                    d_B = -L_bfgs_B(L_BFGS_B.rho, gradient_B_0, L_BFGS_B.S, L_BFGS_B.Y, gamma, Opts);
                else
                    flag_1 = 0;
                    d_B = -gradient_B_0;
                end
                flag_2 = 1;
            else
                flag_1 = 0;
                x1 = x_temp_0;
                ofun_val_1 = ofun_val_0;
                current_Model = convert_vector_to_Model_ROM(Opts.para, x1);
                if i_inner <= 2
                    gradient_diff_B = [];
                    x_diff_B = [];
                end
                break;
            end
        end
    end
end
step_size = alpha;
if i_inner == max_iter
    x1 = x_temp_1;
    current_Model = convert_vector_to_Model_ROM(Opts.para, x1);
    flag_1 = 2;
end
end