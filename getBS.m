function bodysizes=getBS(nof_bins)
    % Using the tree reported in Lee et al. 2014
    % Body size values taken from the tree in Fig 1, 
    % the dates are from the main text and via digitizing Fig 2B
    totalyears=24e7;
    bs_change=[220.7 219 217.4 162.8 46.5 35.8 27 19.1 13.3 9.8 5.5 3.3 1.9 1.2 0.9 0.8 0.8];
    y_change=1e6*[240 235 224 198 174 173.5 173 172 171 170 168.5 167.5 166 165 164 163 0];    
    y_change2=totalyears-y_change;
    times_in_year=linspace(0,totalyears,nof_bins);
    bodysizes=interp1(y_change2,bs_change, times_in_year,'linear');
end