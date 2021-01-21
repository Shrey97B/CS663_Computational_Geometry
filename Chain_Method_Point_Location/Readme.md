Name: Shrey Bhatt

The chain method is used to determine the location of a point in a planar graph structure. Here, it is assumed that the given graph is regular.

The file main.c contains the code in which the chain method is implemented
The files Data1.txt, Data2.txt and Data3.txt contains sample inputs that are provided in console.
The files Result1.txt, Result2.txt and Result3.txt contains outputs of console for the above inputs.
The files Image1.png, Image2.png and Image3.png show the representation of graphs formed in the above inputs respectively.

The input format for the problem is:
N
x1 y1
x2 y2
.
.
.
xN yN
M
p11 p12
p21 p22
p31 p32
.
.
pM1 pM2
Q
q1x q1y
q2x q2y
q3x q3y
.
.
qQx qQy

Explaination: 
Here, N is the number of vertices, after which there are N lines consisting of xi and yi to depict the coordinate of the vertices.
M is the number of edges after which there are M lines consisting of pi1 and pi2.
It indicates that there is an edge between the points having index pi1 and index pi2, which is a zero based indexing to the list of N points listed before.
These indexes provided for edges are based on the original order and not the order after sorting, which is applied in the algorithm.
Q is the number of queries
Each query q gas qix and qiy indicating the coordinates of the query point.

Output Format:
Query point is left to the entire planar subdivision
or
Query point is right to the entire planar subdivision
or
Query point is below the entire planar subdivision
or
Query point is above the entire planar subdivision
or
Query point is just right of the segment ((x1,y1),(x2,y2))
or
Query point is just left of the segment ((x1,y1),(x2,y2))
or
Query point is just on the segment ((x1,y1),(x2,y2))

Explaination: The statements are self explanatory and based on how the super tree nodes are accessed.
