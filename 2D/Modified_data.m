function train_data = Modified_data(train_data, Opts)

train_data.Yu_modified = Opts.para.lambda * train_data.Yu;

end


