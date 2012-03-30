function [new_feature] = make_featurevector_symmetric(featurelist_old)


featurelist_old = swap_this(featurelist_old,1,1,0);
featurelist_old = swap_this(featurelist_old,2,2,1);
featurelist_old = swap_this(featurelist_old,3,3,0);
featurelist_old = swap_this(featurelist_old,4,4,1);
featurelist_old = swap_this(featurelist_old,5,17,0);
featurelist_old = swap_this(featurelist_old,6,18,1);
featurelist_old = swap_this(featurelist_old,7,15,0);
featurelist_old = swap_this(featurelist_old,8,16,1);
featurelist_old = swap_this(featurelist_old,9,13,0);
featurelist_old = swap_this(featurelist_old,10,14,1);
featurelist_old = swap_this(featurelist_old,11,11,0);
featurelist_old = swap_this(featurelist_old,12,12,1);



new_feature = featurelist_old;




end


function [a] = swap_this(a_old,from,to,spiegel)


a_old = a_old';

from_1 = (from-1)*8+1;
from_2 = from_1+7;

to_1 = (to-1)*8+1;
to_2 = to_1+7;

vect_1 = a_old(1,from_1:from_2);
vect_2 = a_old(1,to_1:to_2);

if spiegel == 1
    vect_1=fliplr(vect_1);
    vect_2=fliplr(vect_2);
end

a = a_old;

a(1,from_1:from_2) = vect_2 ;
a(1,to_1:to_2) = vect_1 ;

a = a';



end