% generate binary file for BBFMM3D
Ns = 1e4;
Nf = 1e4;
L = 1;
dof_s = 9;
dof_f = 6;
m = 1;
source = (rand(Ns,3) - 0.5) .* L;
field = (rand(Nf,3) - 0.5) .* L;
q = rand(Ns * dof_s * m,1);
source_binary = reshape(source,Ns * 3,1);
field_binary = reshape(field,Nf*3,1);

fid = fopen('source1e4_test.bin', 'w');
fwrite(fid, source_binary, 'double');
fclose(fid);

fid = fopen('field1e4_test.bin', 'w');
fwrite(fid, field_binary, 'double');
fclose(fid);

fid = fopen('charge1e4_test.bin', 'w');
fwrite(fid, q, 'double');
fclose(fid);
