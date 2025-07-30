function Model = convert_vector_to_Model_ROM(para, x)

x = mat2cell(x, para.partition, 1);

A = reshape(x{1}, para.size_A(2), para.num.A);
A = [A; zeros(para.del, para.num.A)];
A = mat2cell(A, para.size_A(2) + para.del, ones(1,para.num.A))';
[rows, cols] = find(tril(true(para.size_A(1))));
for i = 1:para.num.A
    A{i} = sparse(rows, cols, A{i}, para.size_A(1), para.size_A(1));
    A{i} = tril(A{i}) + tril(A{i}, -1).';
end

temp_B = reshape(x{2}, para.r, para.num.B*para.n_f);
B = mat2cell(temp_B, para.r, para.n_f*ones(1,para.num.B))';

temp_C = reshape(x{3}, para.size_C, para.num.C);
temp_C = mat2cell(temp_C, para.size_C, ones(1,para.num.C))';
row = [para.n_ou*para.rus, para.n_op*para.rp];
C = cell(size(temp_C,1),2);
for j = 1:para.num.C
    temp = mat2cell(full(temp_C{j}), row, 1);
    C{j,1} = reshape(temp{1}, para.n_ou, para.rus);
    C{j,2} = reshape(temp{2}, para.n_op, para.rp);
end

Model.A = A;
Model.B = B;
Model.C = C;

end