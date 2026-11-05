function s = ct_json(name)
%CT_JSON  Read and decode a JSON file from the data/ folder (cached).
persistent cache
if isempty(cache), cache = struct('k', {{}}, 'v', {{}}); end
i = find(strcmp(cache.k, name), 1);
if ~isempty(i), s = cache.v{i}; return; end
fid = fopen(fullfile(ct_root(), 'data', name), 'r', 'n', 'UTF-8');
if fid < 0, error('ct_json:open', 'cannot open data/%s', name); end
txt = fread(fid, inf, '*char')';
fclose(fid);
s = jsondecode(txt);
cache.k{end + 1} = name;
cache.v{end + 1} = s;
end
