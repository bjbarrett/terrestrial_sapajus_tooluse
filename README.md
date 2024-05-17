# terrestrial_sapajus_tooluse
This is a description of the data files in this repo.
## Questions

1. Do monkeys use tools more frequently in areas where they go more to the ground outside tool use context?
  - terrest (out of context)
  - palm density
  - stone availability
  - rio? but not direct, for prediction not causal
  
2. What are the causes of terrestriality in capuchins?
 - tree density (overall trees)
 - distance to river
 - stone tool use


## `nutcracking_data_cropped`
- "timestamp": timestamp of focal with yyyy/mm/dd hh:mm:ss
- "individual": individial names
- "location.l"  : latitude             
- "location_1" : longitude             
- "geometry": spatial info coordinates of grid // 4 corners or centroid                
"corresponding_cell_id_2": cell id number in 22 m grid

## `heights_data`
- "sample_id": `yyyy-mm-dd-_followID_MVIFile_sampleMinute`  i.e. `2022-08-02_0831_MVI_2779_1`
- "focal": `yyyy-mm-dd-_followID_MVIFile` `2022-08-02_0831_MVI_2779`
- "minute": minute in follow where height was recorded, 0-5   
- "timestamp": `yyyy-mm-dd hh:mm:ss` is follow end time i.e `2022-08-02 08:15:00`
- "individual": individual name  
- "sex": male or female         
- "age": adult or juvenile         
- "height": height from floor in meters, 0 on up
- "nutcracking" : 1 yes during follow, 0 no
