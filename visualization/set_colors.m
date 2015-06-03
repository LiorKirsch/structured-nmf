function [clrs,altclrs] = set_colors(shift)
%
%
%
    if(exist('shift')==0), shift = 0;end  
    clrs{1}      =[0.5 0.5 0.99];
    clrs{2}      =[0.31 0.4 0.58];
    clrs{3-shift}=[0.9 0.3 0.0];
    clrs{4-shift}=[1.0 0.7 0.3];
    clrs{5-shift}=[0.2 0.7 0.4];
    clrs{6-shift}=[0.6 0.1 0.05];
    clrs{7-shift}=[0.0 0.5 0.5];
    clrs{8-shift}=[0.5 0.2 0.5];
    clrs{9-shift}=[0.5 0.5 0.2];
    clrs{10-shift}=[0.2 0.5 0.5];
    clrs{11-shift}=[0.9 0.25 0.0];
    clrs{12-shift}=[0.3 0.25 0.9];
    clrs{13-shift}=[0.1 0.9  0.0];
    clrs{14-shift}=[0.9 0.25 0.0];
    clrs{15-shift}=[0.1 0.25 0.9];
    clrs{16-shift}=[0.2 0.7 0.0];
    clrs{17-shift}=[0.4 0.4 0.8];
    
    altclrs{1}  =  [ 0.9501  0.8154  0.0579 ];
    altclrs{2}  =  [ 0.4860  0.7382  0.9099 ]; 
    altclrs{3}  =  [ 0.9311  0.4019  0.1929 ];
    altclrs{4}  =  [ 0.6501  0.5154  0.5579 ];
    altclrs{5}  =  [ 0.4860  0.7382  0.9099 ]; 
    altclrs{6}  =  [ 0.9311  0.4019  0.1929 ];
    altclrs{7}  =  [ 0.4565  0.1355  0.6987 ];
    altclrs{8}  =  [ 0.0185  0.9169  0.6038 ];
    altclrs{9}  =  [ 0.8214  0.4103  0.2722 ];
    altclrs{10} =  [ 0.4447  0.8936  0.1988 ];
    altclrs{11} =  [ 0.4966  0.7271  0.7948 ];
    altclrs{12} =  [ 0.8998  0.3093  0.9568 ];
    altclrs{13} =  [ 0.8216  0.8385  0.5226 ];
    altclrs{14} =  [ 0.6449  0.5681  0.8801 ];
    altclrs{15} =  [ 0.8180  0.3704  0.1730 ];
    altclrs{16} =  [ 0.6602  0.7027  0.9797 ];
    altclrs{17} =  [ 0.3420  0.5466  0.2714 ];
    altclrs{18} =  [ 0.2897  0.4449  0.2523 ];
    altclrs{19} =  [ 0.3412  0.6946  0.8757 ];
    altclrs{20} =  [ 0.5341  0.6213  0.7373 ];
    
    for(i=1:20)
        clrs{17-shift+i} = altclrs{i};
    end
    for(i=21:40)
        clrs{16-shift+i} = altclrs{i-20};
    end
    for(i=41:60)
        clrs{16-shift+i} = altclrs{i-40};
    end