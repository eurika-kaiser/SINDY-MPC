function fix_dottedline(filename)
% Fix the way that the default dotted line and dash-dot lines look terrible
% when exported as eps

fid = fopen(filename,'r');
tempfile = tempname;
outfid = fopen(tempfile,'w');

repeat = 1;
while repeat==1
    thisLine = fgetl(fid);
    iStart = strfind(thisLine,'/DO { [.5');
    if iStart
        thisLine(iStart+7:iStart+8) = '03';
    end
    iStart = strfind(thisLine,'/DD { [.5');
    if iStart
        thisLine(iStart+7:iStart+9) = '1.5';
        thisLine(iStart+10:end+1) = [' ' thisLine(iStart+10:end)];
    end

    if ~ischar(thisLine)
        repeat = 0; 
    else
        fprintf(outfid,'%s\n',thisLine);
    end
end
    
fclose(fid);
fclose(outfid);
copyfile(tempfile, filename);
delete(tempfile);