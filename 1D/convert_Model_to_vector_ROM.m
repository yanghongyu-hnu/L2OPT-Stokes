function result = convert_Model_to_vector_ROM(ROM_Model, para)

result = [];

for j = 1:para.num.A
    temp = ROM_Model.A{j}(tril(true(size(ROM_Model.A{j}))));
    result = [result; temp(1:end-para.del)];
end

for j = 1:para.num.B
    result = [result; ROM_Model.B{j}(:)];
end

for j = 1:para.num.C
    result = [result; ROM_Model.C{j,1}(:)];
    result = [result; ROM_Model.C{j,2}(:)];
end
