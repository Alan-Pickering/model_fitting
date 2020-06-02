function [filename, pathname, filterindex ] = myfilesel(default_file, startpath, in_or_out)
%this just selects a file

cd(startpath);

filterindex=0;
while filterindex==0
    
    if strcmp(in_or_out,'save')
        [filename, pathname, filterindex] = uiputfile( '*.xlsx',  'Select File for Saving Model Fits', default_file);
    elseif strcmp(in_or_out,'read')
        [filename, pathname, filterindex] = uigetfile( '*.xlsx',  'Select File for Reading Task Data', default_file);
    end;
    
    if filterindex==0
        menu('You did not select a filename','ok');
        filterindex=2;
    end; 

end;

end

