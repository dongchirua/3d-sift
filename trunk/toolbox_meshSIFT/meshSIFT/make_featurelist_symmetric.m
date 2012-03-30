function featurelistsym = make_featurelist_symmetric(featurelist)
% featurelistsym = make_featurelist_symmetric(featurelist)
featurelistsym = featurelist;

for i = 1:length(featurelist)
    featurelistsym{1,i}.feature = make_featurevector_symmetric(featurelist{1,i}.feature);
    
end