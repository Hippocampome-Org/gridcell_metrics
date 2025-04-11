% created by Nate Sutton 2025
FID = fopen('test.dat', 'w');
STR =':Frame        1       1 "2002/01/01 07:00:00.000';
fprintf(FID, '%50s \n', STR);
fclose(FID);