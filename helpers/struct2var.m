function struct2var(s)
cellfun(@(n,v) assignin('base',n,v),fieldnames(s),struct2cell(s));
end