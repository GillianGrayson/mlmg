function fn = get_up_figures_path()
if strcmp(getenv('computername'), 'MSI') 
    fn = 'C:/Users/user/Google Drive/mlmg/figures'; 
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS') 
    fn = 'D:/Aaron/Bio/mlmg/figures'; 
else 
    fn = 'C:/Users/user/Google Drive/mlmg/figures'; 
end
end