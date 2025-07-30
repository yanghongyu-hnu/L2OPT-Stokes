function current_Model = recover_from(current_Model, para)

lambda = para.lambda^(1/2);

for i = 1:para.num.C
    current_Model.C{i,1} = current_Model.C{i,1}/lambda;
    current_Model.C{i,2} = current_Model.C{i,2}*lambda;
end

for i = 1:para.num.A
    current_Model.A{i} = current_Model.A{i}*lambda;
end
