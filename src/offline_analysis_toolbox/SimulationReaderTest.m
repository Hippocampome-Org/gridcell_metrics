% created by Nate Sutton 2025
fid = fopen('test.dat', 'r', 'l');
sign = fread(fid, 1, 'int32');
fclose(fid);