function minval = minfunCO2PP(co20s, PPis, DC_Rs, DC_Mus, Fit2datB, m_DIC, m_ph, m_CO2, size_frac_index, amount_calcite_data, d13C_data)  

    % size_fsVE is Matrix F from the manuscript
    size_fsVE = forward_13C(co20s, PPis, DC_Rs, DC_Mus,  Fit2datB, m_DIC, m_ph, m_CO2, size_frac_index, amount_calcite_data);
    minval      = nanmean(nanmean((d13C_data - size_fsVE).^2)); 
    minval = minval + 10^100 * any(PPis <= 0);

end
