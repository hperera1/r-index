### CSA random access data structure:

### what has been added
ri_rasa.hpp: handles CSA data structure construction and querying the data structure.

ri_rasa_tree.hpp: handles tree construction and traversal at query time.

(within r_index.hpp) CSA construction happens after building phi. In order to query the r-index an additional function called query_csa has been added which takes in an index (sa_i) and the SA samples needed (ssa).

### what hasn't been tested 

ri_rasa.hpp:
build_rads_phi_inv() has not been tested as I used Phi while writing the code. 

find_cycles() hasn't been properly tested but it hasn't caused any problems yet.

### how to test query results 

ri-query & rasa-query have useful help() functions that provide the best information I think on how to use them.

neither of them output anything unless you make them as well, just thought I'd let you know.

ri-buildfasta builds the CSA after building the phi function. The flags -b bigbwt were always used.

ri-query, when used to query up to a particular value, will output those values to a text file in the /rasa_tests folder supposing lines around 185 are uncommented. 

rasa-query can then be used to scan this file (line 70) and the samples against the samples that it retrieves as it iterates and queries each sample.

If rasa-query retrieves an incorrect sample you can see which sample it was looking for, which was retrieved, and at which index. I would use the index and set a break point when 'i' was the index that was giving problems.
