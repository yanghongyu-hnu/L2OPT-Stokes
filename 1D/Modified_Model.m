function result = Modified_Model(Model, Opts)

lambda = Opts.para.lambda^(1/2);

[A, B, C, num] = read_Model(Model);

result.A = cell(Opts.para.num.A,1);

for i = 1:num.A
    result.A{i} = A{i}/lambda;
end

result.B = Model.B;

result.C = cell(Opts.para.num.C,Opts.para.num_output);
for i = 1:Opts.para.num.C
    result.C{i,1} = Model.C{i,1}*lambda;
    result.C{i,2} = Model.C{i,2}/lambda;
end

end


