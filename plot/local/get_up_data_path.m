function fn = get_up_data_path()
if strcmp(getenv('computername'), 'MSI') 
    fn = 'D:/YandexDisk/Work/mlmg'; 
elseif strcmp(getenv('computername'), 'DESKTOP-4BEQ7MS') 
    fn = 'C:/Users/User/YandexDisk/mlmg'; 
else 
    fn = 'E:/YandexDisk/Work/mlmg'; 
end 
end