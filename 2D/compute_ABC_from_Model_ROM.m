function [A_hat, B_hat, C_hat] = compute_ABC_from_Model_ROM(ROM_Model, para, ROM_coe_fun)

A = ROM_Model.A; B = ROM_Model.B; C = ROM_Model.C;

A_hat = 0; B_hat = 0; C_hat = cell(1,para.num_output);

for k = 1:para.num_output
    C_hat{1,k} = 0;
end

for j = 1:para.num.A
    A_hat = sparse(A_hat + sparse(feval(ROM_coe_fun.A{j}, para.int_para) * A{j}));
end

for j = 1:para.num.B
    B_hat = sparse(B_hat + sparse(feval(ROM_coe_fun.B{j}, para.int_para) * B{j}));
end

for j = 1:para.num.C
    for k = 1:para.num_output
        C_hat{1,k} = sparse(C_hat{1,k} + sparse(feval(ROM_coe_fun.C{j}, para.int_para) * C{j,k}));
    end
end


end