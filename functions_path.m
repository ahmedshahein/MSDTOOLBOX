warning off all;
quantizer

clear
clc
close all

fprintf('#####################\n')
fprintf('Welcome to MSDTOOLBOX\n')
fprintf('#####################\n')

fprintf('\n')
fprintf('Ahmed Shahein\n')
fprintf('Fritz Huettinger Chair of Microelectronics\n')
fprintf('Department of Microsystems Engineering - IMTEK\n')
fprintf('University of Freiburg\n')
fprintf('Georges-Koehler-Allee 102\n')
fprintf('79110 Freiburg\n')
fprintf('Germany\n')


fprintf('\n')
if (ispc),
    today = date;
    fprintf('DATE ... %s \n',today);
end
if isunix,
    ! date '+DATE: %d.%m.%y%nTIME: %H:%M:%S'
end
fprintf('\n')

dir = pwd;
path(path, dir)
cd ./aFunctions
dir = pwd;
path(path, dir)
cd ../bFunctions
dir = pwd;
path(path, dir)
cd ../iPatterns
dir = pwd;
path(path, dir)
cd ../

chkFolder = exist('eCoefficients','dir');
if chkFolder == 7,
    cd eCoefficients
    dir = pwd;
    path(path, dir)
    cd ../
else
    mkdir('eCoefficients')
    cd eCoefficients
    dir = pwd;
    path(path, dir)
    cd ../
end

fprintf('###########################################\n')
fprintf('CONGRATULATIONS: MSDTOOLBOX is ready to use\n')
fprintf('###########################################\n')
fprintf('\n')

fprintf('It is recommended to read the README\n')
fprintf('\n')
replyREADME = input('Do you want to open the README? Y/N [Y]: ', 's');
if strcmp(replyREADME, 'Y'), 
    open README
else
    replyREADME = 'N';
end

fprintf('\n')
reply = input('Do you want start the HELP? Y/N [Y]: ', 's');
if isempty(reply) || strcmp(reply, 'Y'), 
    cd doc
    web index.html
    cd ../
else
    reply = 'N';
end