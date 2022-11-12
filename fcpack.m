function fcpack(dst)
%FCPACK     Pack fcmclt release into a zip file for delivery.
%           fcpack(dst) generates the necessary mex files, packs the
%           fcmclt release into a ZIP file and places it inside the
%           destination directory.
%
%           The dst is a string containing the directory to store the ZIP file.
%
%           The destination directory must be an absolute path
%           and must not already exist since it will be created
%           by this script.
%
%           The documentation and license file will be packed into a zip file as well.

% MIT License
%
% Copyright (c) 1995-2022 Bogdan Kosanovic
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Checking if destination is missing '\' or '/' at the end
ind1 = findstr('\',dst);
ind2 = findstr('/',dst);
if length(ind1) < 1, ind1 = 0; end
if length(ind2) < 1, ind2 = 0; end
if length(dst) ~= max(ind1) && length(dst) ~= max(ind2);
  dst = [dst,'\'];
end

% Find the location of this script
source = mfilename('fullpath');
ind = find(source == '\');
ind = max(ind);
source = source(1:ind);
source = lower(source);

% Check if destination exists
dst = lower(dst);
if isdir(dst);
  error(sprintf('Destination directory %s already exists. Provide different directory',dst));
end

cd([source, 'src']);
[dumm,VERSION_FCMCLT] = fcver();
cd(source);

FILE_NAME = ['fcmclt_', VERSION_FCMCLT, '_w64'];
ZIP_FILE = [FILE_NAME, '.zip'];

disp('... creating destination folder'); drawnow;
% Change current directory to build files 
oldcurdir = cd;             % current directory saved
csrc = [source, 'src\'];    % this is where the C files are
mkdir([dst,'temp']);

disp('... copying files to destination'); drawnow;

% Copy source code, etc. to destination (so that ZIP would work properly!)
copyfile([source,'*.md'],[dst, 'temp']);
copyfile([source,'src\*.m'],[dst, 'temp\src']);
copyfile([source,'src\*.md'],[dst, 'temp\src']);
copyfile([source,'src\*.c'],[dst, 'temp\src']);

% Go to temp src folder to build mex files
cd([dst,'temp\src']);

% Build mex-files
disp('... building fcmclt mex files'); drawnow;
mex extfpm.c
mex fcmc.c
mex gkfcmc.c

% Go back to temp folder root
cd([dst, 'temp']);

% Zip the files
fn=1;     files{fn,:} = '*';

disp('... packing files into the archive'); drawnow;
zip(ZIP_FILE,files);
cd(dst);                                    % go back up to desintation folder
copyfile(['temp\', ZIP_FILE],ZIP_FILE);

% Delete temporary files
rmdir('temp','s');

% Go back to where you started from
cd(oldcurdir);

disp('Done.');

% nothing past this point

