function maindir=SetBEMPaths
% function maindir=SetBEMPaths
% sets Matlab paths for the BEM library

%type here the path, where the library is located
maindir='/Volumes/grillo/Benutzer/dp032/Tmp/Matti/BEMlibrary_100930';

paths=genpath(maindir);
addpath(paths);
