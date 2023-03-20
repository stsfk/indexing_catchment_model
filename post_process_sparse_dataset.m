load("sparse_out.mat")


sparse_index = eval_grid(:,[1 3]);
sparse_parameter = eval_grid(:,2);
sparse_R = eval_grid(:,4);

writecell(sparse_index,'sparse_index.csv')
writecell(sparse_parameter,'sparse_parameter.csv');
writecell(sparse_R,'sparse_R.csv');