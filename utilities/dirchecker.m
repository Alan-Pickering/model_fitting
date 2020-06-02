function dirchecker(currDir,homeDirectory)

%checks the directory a script is launched from
if strcmp(currDir,homeDirectory)
  %all is good
else
  disp('You were in the wrong directory, as below:-')
  disp(currDir);
  disp('Restoring to the correct home directory.');
  cd(homeDirectory);  
  disp('Hit any key to continue');
  pause;
  clc;
end;



end

