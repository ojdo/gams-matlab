% fuelstation.m using interface class FS
f = FS;
% set path to gams executable manually, if gams is not on your system path
f.path.gams = 'C:/GAMS/win64/24.5/gams.exe';
f.writeInputs;
f.run;

if f.status == 0
    % if run was successful, read result and plot variable x 
    f.readResults;
    f.plot;
else
    % if not, display error message
    display(char(f.result{:}))
end
