3yrprj
======

How is a set represented?


There will be two definitions, an asset set, where you define your assets as an array(could be removed in future). 
And a member-set/subset where you define which members this set contains from set A.
A member-set/subset will be defined by n numbers of intigers where one integer represent 32 members of that set.

If the set A(represented as an array) =  apple, mapple, banana, potato where apple have the index of 0, mapple 1... in the array.
A subset/member-set of A containing banana and apple, the integer representing which members contained would have a value of 5 as 0101 in binary is equal to 5.
So the least significant digit of the member-set represent the first element in the array A.

The reason behind this design choise is that 

Speed Graph

 ![Alt text](https://docs.google.com/spreadsheet/oimg?key=0Au1A39JI-BwHdE1YRHFWVlpkXzZKWGxNS0tHMmRmVHc&oid=0&zx=93h5htrgmj8j "Speed graph")