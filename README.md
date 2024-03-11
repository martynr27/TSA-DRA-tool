# TSA-DRA-tool
The Temporary Storage Area Drainage Rate Analysis (TSA-DRA) tool is a data-based mechanistic method that can assess the functioning of nature-based solutions (NbS) and only requires rainfall and water level / volume data.

The tool provides a systematic approach for characterising the functioning of a wide range of temporary storage area (TSA) types and sizes. The tool can also explore time-variable TSA functioning and provides useful soil infiltration estimates for modelling.

Method:
This script extracts recession periods based on user input and then creates a master recession curve (MRC) to provide a general description for TSA drainage over a longer time period. A segmented linear model is then fitted to the MRC to describe the TSA drainage rate. User input is based on site knowledge and the desired MRC analysis. 

Data requirements:
TSA data | 15 min | Headings: date_time, temp, depth 
Precipitation | hourly | Headings: date_time, mm_hr
