addpath('/Volumes/bordeaux/IBT10.10/matlab/common/vtkToolbox/MATLAB:');

merge = vtkRead('/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/vtk_files/tissue_ints_60_f1.vtk');

writeCARP(merge,'/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/8x8/tissue_ints_60_f1/tissue_ints_60_f1');



function writeCARP(inStruct,filename)

 

pts_size = size(inStruct.points,1);

cell_size = size(inStruct.cells,1);

 

fileID = fopen([filename '.pts'],'w');

fprintf(fileID,'%i\n',pts_size);

fprintf(fileID,'%f %f %f\n',inStruct.points.');

fclose(fileID);

 

if inStruct.cellTypes(1) == 5

    

    fileID = fopen([filename '.elem'],'w');

    fprintf(fileID,'%i\n',cell_size);

    fprintf(fileID,'Tr %i %i %i %i\n',[inStruct.cells.'-1 inStruct.cellData.elemTag.']);

    fclose(fileID);

    

    if size(inStruct.cellData.fiber,2) == 3

    elseif size(inStruct.cellData.fiber,2) == 6

    else

    end

elseif inStruct.cellTypes(1) == 10

    

    fileID = fopen([filename '.elem'],'w');

    fprintf(fileID,'%i\n',cell_size);

    fprintf(fileID,'Tt %i %i %i %i %i\n',[inStruct.cells.'-1; inStruct.cellData.elemTag.']);

    fclose(fileID);

    

    if size(inStruct.cellData.fiber,2) == 3

        

        fileID = fopen([filename '.lon'],'w');

        fprintf(fileID,'1\n');

        fprintf(fileID,'%f %f %f\n',inStruct.cellData.fiber.');

        fclose(fileID);

        

    elseif size(inStruct.cellData.fiber,2) == 6

        fileID = fopen([filename '.lon'],'w');

        fprintf(fileID,'2\n');

        fprintf(fileID,'%f %f %f %f %f %f\n',inStruct.cellData.fiber.');

        fclose(fileID);

    else

        fileID = fopen([filename '.lon'],'w');

        fprintf(fileID,'1\n');

        fprintf(fileID,'%f %f %f\n',zeros(cell_size,3));

        fclose(fileID);

    end

else

    disp('Element type not yet define')

end

 

end
