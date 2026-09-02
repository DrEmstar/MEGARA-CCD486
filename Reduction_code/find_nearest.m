function nearest = find_nearest(vector,value)
test=vector-value;
nearest=find( abs(test) == min(abs(test)) );
if numel(nearest) > 1
    nearest=nearest(floor(numel(nearest)/2));
end
