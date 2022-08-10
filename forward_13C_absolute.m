
function size_fsVE = forward_13C_absolute(co20s, PPis, DC_Rs, DC_Mus, Fit2datB, m_DIC, m_ph, m_CO2, size_frac_index, amount_calcite_data)
    fPPs = repmat(PPis, [length(co20s), 1]);

    DIC1 = interp1(m_CO2, m_DIC, co20s);     % interpolate 
    ph1 = interp1(m_CO2, m_ph, co20s); 
    fDIC = repmat(DIC1, [1, length(PPis)]);
    fph = repmat(ph1, [1, length(PPis)]);
    
    d13Cmod     =  ISOMOD_analytical_fast(fDIC(:), fph(:), 0, ...            % env
                                    DC_Rs(:), DC_Mus(:), fPPs(:), 1, ...       % bio
                                    Fit2datB);  %                        % consts
    dCfit  =  NaN(size(fDIC)) ;
    dCfit(:)    =   d13Cmod(:,4);
    size_fs = (size_frac_index * (amount_calcite_data .* dCfit'))  ./  (size_frac_index * (amount_calcite_data));
    size_fsVE = size_fs;
    
    
end

