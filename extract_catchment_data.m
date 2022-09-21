function [catchment_data] = extract_catchment_data(catchment_datas,catchment_id)
%   extract_catchment_data: extract catchment data from catchment_datas based
%on catchment_id
catchment_data = catchment_datas{catchment_id,1};
end

