Name: Shrey Bhatt

The given program implements the regularization of planar graph using plane sweep method.

The file main.c contains the code in which the regularization method is implemented
The files Input_01.txt, Input_02.txt, Input_03.txt, Input_04.txt contains sample inputs that are provided in console.
The files Output_01.txt, Output_02.txt, Output_03.txt, Output_04.txt contains outputs of console for the above inputs.
The files Reg_01.png, Reg_02.png, Reg_03.png, Reg_04.png show the representation of graphs formed in the above inputs respectively.

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

Explaination: 
Here, N is the number of vertices, after which there are N lines consisting of two values.
xi and yi depict the coordinate of the vertices.
M is the number of edges after which there are M lines consisting of pi1 and pi2.
It indicates that there is an edge between the points having index pi1 and index pi2, which is a one based indexing to the list of N points listed before.
These indexes provided for edges are based on the original order and not the order after sorting, which is applied in the algorithm.

Output Format:

No new edges have been added

or

New edges added between points:

(x11,y11) (x12,y12)

(x21,y21) (x22,y22)

(x31,y31) (x32,y32)

.
.
.

Explaination: If there are no new edges added in both the pass, the output is "No new edges have been added".
If edges are added, the other output is displayed.
Each line indicates an edge added and contains pair of points.
Each pair of points (xi1,yi1) (xi2,yi2) indicates the co-ordinates of points between which the edges are added.
i.e. An edge is added between points with the coordinates (xi1,yi1) (xi2,yi2).
