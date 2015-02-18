function idx = search(x, data)

    if ~isscalar(x)
        idx = arrayfun(@(y)search(y,data),x);
        return
    end

    if (length(data) < 1e5)
        if (x >= data(end))
            idx = length(data);
        else
            idx = find(x < data, 1, 'first');
        end
        return;
    end

    imin = 1;
    imax = length(data);
    
    while (imin <= imax)
        idx = ceil((imin + imax)/2);
        if (x <= data(idx))
            imax = idx - 1;
        else
            imin = idx + 1;
        end
    end
    idx = imin;

end