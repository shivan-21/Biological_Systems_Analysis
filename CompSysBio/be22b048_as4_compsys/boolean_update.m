function [new_state, cell_fate] = boolean_update(current_state)
    % Unpack current states
    P53 = current_state(1);
    MDM2 = current_state(2);
    MYC = current_state(3);
    RB = current_state(4);
    
    % Update rules (Part A)
    new_P53 = MYC && ~MDM2;  % P53 activated by MYC
    new_MDM2 = MDM2; % MDM2 activated by itself
    new_MYC = ~RB;  % MYC inhibited by RB
    new_RB = ~MYC;  % RB inactivated by MYC
    
    % Determine cell fate (Part D)
    apoptosis = P53 && ~MDM2;
    proliferation = MYC && ~RB;
    
    new_state = [new_P53, new_MDM2, new_MYC, new_RB];
    cell_fate = '';
    if apoptosis, cell_fate = [cell_fate 'A']; end
    if proliferation, cell_fate = [cell_fate 'P']; end
    if isempty(cell_fate), cell_fate = 'None'; end
end
