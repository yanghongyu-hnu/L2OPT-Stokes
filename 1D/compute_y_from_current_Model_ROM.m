function result = compute_y_from_current_Model_ROM(ROM_Model, para, ROM_coe_fun, Xi)

N_xi = size(Xi, 1);

result.Yu_hat = sparse(N_xi, para.n_ou * para.n_f);
result.Yp_hat = sparse(N_xi, para.n_op * para.n_f);

for i = 1:N_xi
    para.int_para.xi = Xi(i);
    [A_hat, B_hat, C_hat] = compute_ABC_from_Model_ROM(ROM_Model, para, ROM_coe_fun);
    x_hat = A_hat\B_hat;

    u_hat = x_hat(1:para.rus,:);
    p_hat = x_hat(para.rus+1:end,:);

    yu_hat = C_hat{1} * u_hat;
    yp_hat = C_hat{2} * p_hat;

    result.Yu_hat(i,:) = yu_hat(:)';
    result.Yp_hat(i,:) = yp_hat(:)';
end
end
