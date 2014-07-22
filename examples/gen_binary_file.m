% Example of generating binary file for BBFMM3D

Ns = 1e4; % Number of sources in simulation cell
Nf = 1e4; % Number of field points in simulation cell
L = 1;    % Length of simulation cell (assumed to be a cube)
dof_s = 9;% Dimension of the tensor kernel
dof_f = 6;
m = 1;    % Number of columns of the vector
source = (rand(Ns,3) - 0.5) .* L; % Positions of source points
field = (rand(Nf,3) - 0.5) .* L;  % Positions of field points
q = rand(Ns * dof_s,m);       % charge

source_binary = reshape(source,Ns * 3,1);
field_binary = reshape(field,Nf * 3,1);
q_binary = reshape(q,Ns * dof_s * m,1);
fid = fopen('source1e4_test.bin', 'w');
fwrite(fid, source_binary, 'double');
fclose(fid);

fid = fopen('field1e4_test.bin', 'w');
fwrite(fid, field_binary, 'double');
fclose(fid);

fid = fopen('charge1e4_test.bin', 'w');
fwrite(fid, q, 'double');
fclose(fid);
