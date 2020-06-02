function [filename, pathname, filterindex ] = myfilesel_cs(default_file, startpath, in_or_out)
%this just selects a file

cd(startpath);

filterindex=0;
while filterindex==0
    
    if strcmp(in_or_out,'save')
        [filename, pathname, filterindex] = uiputfile( '*.xlsx',  'Select File for Saving in Excel', default_file);
    elseif strcmp(in_or_out,'read')
        [filename, pathname, filterindex] = uigetfile( '*.mat',  'Select File for Reading Task Data Mat files', default_file);
    end;
    
    if filterindex==0
        menu('You did not select a filename','ok');
        filterindex=2;
    end; 

end;

end

