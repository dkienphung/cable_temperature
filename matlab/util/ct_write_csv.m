function ct_write_csv(path, header, rows)
%CT_WRITE_CSV  Write a cell array (rows x columns; numbers or char) with a header line.
fid = fopen(path, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
for i = 1:size(rows, 1)
    parts = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{i, j};
        if ischar(v)
            parts{j} = ['"' strrep(v, '"', '""') '"'];
        elseif islogical(v)
            parts{j} = sprintf('%d', v);
        elseif isempty(v) || (isnumeric(v) && isnan(v))
            parts{j} = '';
        else
            parts{j} = sprintf('%.6g', v);
        end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
fclose(fid);
end
