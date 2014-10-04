[pathstrmain, name, ext] = fileparts(mfilename('fullpath'));
currentfolder = pwd;
mayrfolder = ([pathstrmain,'/K-operator-mayr']);
cd(mayrfolder);
try
    % It seems that the following does the job from within matlab.
    % For a more optimized computation you should compile manually with -03
    % like depicted below.
    mex -I"./eigen3" CFLAGS='$CFLAGS -Wall -fPIC' mx_bem_dlp.c bem_common.cpp bem_matlab.cpp -output buildK
catch
    cd(currentfolder);
    error('The compilation from within MATLAB failed, consult this file for more details.');
%% IF THE ABOVE FAILS:
% The following is a bash script we did use before
% You have to enter the location of your mex-compiler after "MEX="
% Afterwards you can use it as a bash script.
% You can also try to do the g++ compilation from shell and mex compilation
% from MATLAB.

% MEX=/Applications/MATLAB_R2012b.app/bin/mex  # PATH TO MATLAB-MEX-COMPILER
%
% DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"; # DIRECTORY OF THIS FILE
% EIGENPATH=$DIR/eigen3 # PATH TO EIGEN-SOURCECODE - Download eigen3 from http://eigen.tuxfamily.org/
% ## COMPILATION
% cd "$DIR";
% g++ -c -O3 -fPIC -I"$EIGENPATH" bem_common.cpp -o bem_common.o
% g++ -c -O3 -fPIC -I"$EIGENPATH" bem_matlab.cpp -o bem_matlab.o
% "$MEX" -output buildK mx_bem_dlp.c bem_common.o bem_matlab.o
end
cd(currentfolder);
