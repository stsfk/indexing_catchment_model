load("out.mat")

parameters = cell2mat(out(2,3));
writematrix(parameters,'dense_parameters.csv')

R = cell2mat(out(2:end,4));
writematrix(R,'dense_R.csv')