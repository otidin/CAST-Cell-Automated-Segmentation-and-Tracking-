if doAppend
appendvec=[appendvec; fname];
[~,idxy]=unique(appendvec );
appendvec=appendvec(idxy,:)
end
