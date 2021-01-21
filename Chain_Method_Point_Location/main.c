//Author Shrey Bhatt 20111060
#include <stdio.h>
#include <stdlib.h>
#include<math.h>

//Structure representing a point
struct point{
    double xcoord;
    double ycoord;
    int Ival;
    struct edgeListNode* outL;
    struct edgeListNode* inL;
    struct edgeListNode* outLTail;
    struct edgeListNode* inLTail;
};

//Structure representing a edge
struct edge{
    struct point* p1;
    struct point* p2;
    int edgNum;

};

//Structure used to store edge list for incoming and outgoing list
struct edgeListNode{
    struct edge *ed;
    struct edgeListNode* next;
};

//Structure to represent super tree node
struct superTreeNode{
    struct edge **chNode;
    int cindex;
    struct superTreeNode* left;
    struct superTreeNode* right;
    struct superTreeNode* par;
};

//n number of vertices and m number of edges
int n,m;
struct point **vertices;
struct edge* edges;
int numchains;

//comparator function for vertex sorting
int vcomp(const void* v1, const void* v2){
    const struct point **ver1 = (const struct point **) v1;
    const struct point **ver2 = (const struct point **) v2;
    double x1=(*ver1)->xcoord;
    double x2=(*ver2)->xcoord;
    double y1=(*ver1)->ycoord;
    double y2=(*ver2)->ycoord;

    if(y1!=y2){
        if(y1<y2){
		return -1;
	}
	if(y1>y2){
		return 1;
	}
	if(y1==y2){
		return 0; //should not reach here
	}
    }
    else{
        if(x2<x1){
		return -1;
	}
	if(x2>x1){
		return 1;
	}
	if(x1==y2){
		return 0; //should not reach here
	}
    }
}

//ATan2 value for calculating angle
double getATan2Val(struct point* p1,struct point* p2){
    double ydiff = p2->ycoord-p1->ycoord;
    double xdiff = p2->xcoord-p1->xcoord;
    return atan2(ydiff,xdiff);
}


struct edgeListNode* SortedMergeOut(struct edgeListNode* a, struct edgeListNode* b)
{
    struct edgeListNode* result = NULL;

    if (a == NULL)
        return (b);
    else if (b == NULL)
        return (a);

    struct edge* aed = a->ed;
    struct edge* bed = b->ed;

    double atan2v = getATan2Val(aed->p1,aed->p2);
    double btan2v = getATan2Val(bed->p1,bed->p2);
    if (atan2v >= btan2v) {
        result = a;
        result->next = SortedMergeOut(a->next, b);
    }
    else {
        result = b;
        result->next = SortedMergeOut(a, b->next);
    }
    return (result);
}

struct edgeListNode* SortedMergeIn(struct edgeListNode* a, struct edgeListNode* b)
{
    struct edgeListNode* result = NULL;

    if (a == NULL)
        return (b);
    else if (b == NULL)
        return (a);

    struct edge* aed = a->ed;
    struct edge* bed = b->ed;

    double atan2v = getATan2Val(aed->p1,aed->p2);
    double btan2v = getATan2Val(bed->p1,bed->p2);
    if (atan2v <= btan2v) {
        result = a;
        result->next = SortedMergeIn(a->next, b);
    }
    else {
        result = b;
        result->next = SortedMergeIn(a, b->next);
    }
    return (result);
}

void listSplit(struct edgeListNode* head,
                    struct edgeListNode** a, struct edgeListNode** b)
{
    struct edgeListNode* fast;
    struct edgeListNode* slow;
    slow = head;
    fast = head->next;

    while (fast != NULL) {
        fast = fast->next;
        if (fast != NULL) {
            slow = slow->next;
            fast = fast->next;
        }
    }
    *a = head;
    *b = slow->next;
    slow->next = NULL;
}

//function to merge outgoing list
void mergesortoutl(struct edgeListNode** head_ref){
    struct edgeListNode* head = *head_ref;
    struct edgeListNode* a;
    struct edgeListNode* b;

    if ((head == NULL) || (head->next == NULL)) {
        return;
    }

    listSplit(head, &a, &b);

    mergesortoutl(&a);
    mergesortoutl(&b);

    *head_ref = SortedMergeOut(a, b);
}

//function to merge incoming list
void mergesortinl(struct edgeListNode** head_ref){
    struct edgeListNode* head = *head_ref;
    struct edgeListNode* a;
    struct edgeListNode* b;

    if ((head == NULL) || (head->next == NULL)) {
        return;
    }

    listSplit(head, &a, &b);

    mergesortinl(&a);
    mergesortinl(&b);

    *head_ref = SortedMergeIn(a, b);
}

//sorting vertices
void sortVertices(struct point** vs){
    qsort((void*) vs,n,sizeof(vs[0]),vcomp);
}

//construction of super tree
struct superTreeNode* recConstrSuperTree(int l,int r,struct superTreeNode* par){
    if(l>r){
        return NULL;
    }
    int countn = r - l+1;
    int numbit=0;
    while(countn>0){
        countn = countn>>1;
        numbit++;
    }
    int nval = l + (1<<(numbit-1)) - 1;
    //printf("Constructing for %d in %d %d\n",nval,l,r);
    struct superTreeNode* root = (struct superTreeNode*) malloc(sizeof(struct superTreeNode));
    root->cindex = nval;
    root->chNode = NULL;
    root->par = par;

    root->left = recConstrSuperTree(l,nval-1,root);
    root->right = recConstrSuperTree(nval+1,r,root);
    return root;

}

/*
void travSuperTree(struct superTreeNode* root,int *numedges){
    if(root==NULL){
        return;
    }
        //printf("Processing for %d\n",root->cindex);
        int x = numedges[root->cindex];
        for(int i=0;i<x;i++){
            struct edge* ed = root->chNode[i];
            //printf("%d\t",ed->edgNum);
        }
        //printf("\n");
        travSuperTree(root->left,numedges);
        travSuperTree(root->right,numedges);
}
*/

struct superTreeNode* findHeight(struct superTreeNode* root,int val,int* height){
    if(val==root->cindex){
        return root;
    }
    if(val<root->cindex){
        (*height) = (*height) + 1;
        return findHeight(root->left,val,height);
    }
    else{
        (*height) = (*height) + 1;
        return findHeight(root->right,val,height);
    }
}


//finding LCA for two chains in super tree
struct superTreeNode* findLCA(struct superTreeNode* root,int leftv,int rightv){
    //printf("left %d right %d %d %d\n",leftN->cindex,rightN->cindex,hl,hr);
    if(root==NULL){
        return NULL;
    }

    if(root->cindex>leftv && root->cindex>rightv){
        return findLCA(root->left,leftv,rightv);
    }

    if(root->cindex<leftv && root->cindex<rightv){
        return findLCA(root->right,leftv,rightv);
    }

    return root;

}

void initializeChainNodes(struct superTreeNode* root,int *edgInChain){
    if(root==NULL){
        return;
    }

    int cind = root->cindex;
    int numedg = edgInChain[cind];
    root->chNode = (struct edge**) (malloc(numedg * sizeof(struct edge*)));
    initializeChainNodes(root->left,edgInChain);
    initializeChainNodes(root->right,edgInChain);

}

//assigning an edge to a chain node of super tree
void assignToChainNum(struct superTreeNode* root,struct edge* edgptr,int chainNum,int currind){

    if(chainNum==root->cindex){
            root->chNode[currind] = edgptr;
            return;
    }

    if(chainNum<root->cindex){
        assignToChainNum(root->left,edgptr,chainNum,currind);
    }
    else{
        assignToChainNum(root->right,edgptr,chainNum,currind);
    }
}

int main()
{
    int p1v, p2v;
    scanf("%d",&n); //number of vertices
    vertices = (struct point **) malloc(n * sizeof(struct point*));

    //Taking vertices/points input
    for(int i=0;i<n;i++){
        struct point* vertex = (struct point*) malloc(sizeof(struct point));
        scanf("%lf",&(vertex->xcoord));
        scanf("%lf",&(vertex->ycoord));
        vertices[i] = vertex;
    }

    scanf("%d",&m);
    edges = (struct edge *) malloc(m * sizeof(struct edge));

    //Accepting input of edges described by zero-based index of two points/vertex
    for(int i=0;i<m;i++){
        scanf("%d",&p1v);   //0-based index of vertex
        scanf("%d",&p2v);
        edges[i].p1 = vertices[p1v];
        edges[i].p2 = vertices[p2v];
        edges[i].edgNum = i;
    }

    /*
    for(int i=0;i<m;i++){
        printf("Edge %d P1 and P2: %lf %lf ,%lf %lf\n",i,edges[i].p1->xcoord,edges[i].p1->ycoord,edges[i].p2->xcoord,edges[i].p2->ycoord);
    }
    */

    //struct point* p4 = vertices[3];
    //printf("%lf %lf",p4->xcoord,p4->ycoord);
    sortVertices(vertices); //sorting the vertices such that if i<j, yi<yj or if(yi==yj) xi>xj
    //printf("%lf %lf\n",p4->xcoord,p4->ycoord);
    for(int i=0;i<n;i++){
        //printf("%lf %lf\n",vertices[i]->xcoord,vertices[i]->ycoord);
        vertices[i]->Ival = i;
        vertices[i]->outL = NULL;
        vertices[i]->outLTail = NULL;
        vertices[i]->inL = NULL;
        vertices[i]->inLTail = NULL;
    }


    for(int i=0;i<m;i++){

        //keeping p1 as lower index vertex and p2 as higher indexed vertex
        if(edges[i].p1->Ival>edges[i].p2->Ival){
            struct point* temp = edges[i].p1;
            edges[i].p1 = edges[i].p2;
            edges[i].p2 = temp;
        }
        //printf("%lf %lf %lf %lf\n",edges[i].p1->xcoord,edges[i].p1->ycoord,edges[i].p2->xcoord,edges[i].p2->ycoord);

        //inserting the edge to outgoing list of p1
        struct edgeListNode* olt = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
        olt->ed = &edges[i];
        olt->next = NULL;
        if(edges[i].p1->outLTail==NULL){
            edges[i].p1->outL = olt;
            edges[i].p1->outLTail = olt;
        }
        else{
            edges[i].p1->outLTail->next = olt;
            edges[i].p1->outLTail = olt;
        }

        //inserting the edge to incoming list of p2
        struct edgeListNode* ilt = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
        ilt->ed = &edges[i];
        ilt->next = NULL;
        if(edges[i].p2->inLTail==NULL){
            edges[i].p2->inL = ilt;
            edges[i].p2->inLTail = ilt;
        }
        else{
            edges[i].p2->inLTail->next = ilt;
            edges[i].p2->inLTail = ilt;
        }
    }

    for(int i=0;i<n;i++){

        struct edgeListNode* ilh = vertices[i]->inL;
        struct edgeListNode* olh = vertices[i]->outL;
        mergesortinl(&ilh); //sorting the edges in the incoming list to be counter clockwise
        vertices[i]->inL = ilh;
        mergesortoutl(&olh); //sorting the edges in the outgoing list to be counter clockwise
        vertices[i]->outL = olh;
        struct edgeListNode* pr = NULL;
        while(olh!=NULL){
            //printf("%d ",olh->ed->p2->Ival);
            pr = olh;
            olh = olh->next;
        }
        vertices[i]->outLTail = pr;
        //printf("\n");
        pr=NULL;
        while(ilh!=NULL){
            //printf("%d ",ilh->ed->p1->Ival);
            pr=ilh;
            ilh = ilh->next;
        }
        vertices[i]->inLTail = pr;
        //printf("\n");

    }

    //weight balancing procedure
    int weightEdge[m];
    int weightIn[n];
    int weightOut[n];
    int inDeg[n];
    int outDeg[n];

    //initializtion
    for(int i=0;i<m;i++){
        weightEdge[i]=1;
    }

    for(int i=0;i<n;i++){
        inDeg[i] = 0;
        outDeg[i] = 0;
        struct edgeListNode* inh = vertices[i]->inL;
        while(inh!=NULL){
            inDeg[i]++;
            inh = inh->next;
        }
        struct edgeListNode* outh = vertices[i]->outL;
        while(outh!=NULL){
            outDeg[i]++;
            outh = outh->next;
        }
    }

    for(int i=1;i<n-1;i++){
        weightIn[i]=0;
        struct edgeListNode* inh = vertices[i]->inL;
        struct edge* leftOutEdge = vertices[i]->outL->ed;
        while(inh!=NULL){
            weightIn[i] = weightIn[i] + weightEdge[inh->ed->edgNum];
            inh = inh->next;
        }
        if(weightIn[i]>outDeg[i]){
            weightEdge[leftOutEdge->edgNum] = weightIn[i] - outDeg[i] + 1;
        }
    }

    /*
    printf("Weights after first pass\n");
    for(int i=0;i<m;i++){
        printf("%d\n",weightEdge[i]);
    }*/

     for(int i=n-1;i>=1;i--){
        weightOut[i]=0;
        struct edgeListNode* outh = vertices[i]->outL;
        struct edge* leftInEdge = vertices[i]->inL->ed;
        while(outh!=NULL){
            weightOut[i] = weightOut[i] + weightEdge[outh->ed->edgNum];
            outh = outh->next;
        }
        if(weightOut[i]>weightIn[i]){
            weightEdge[leftInEdge->edgNum] = weightOut[i] - weightIn[i] + weightEdge[leftInEdge->edgNum];
        }
    }

    /*
    printf("Weights after second pass\n");
    for(int i=0;i<m;i++){
        printf("%d\n",weightEdge[i]);
    }*/

    int vertL[n];
    int vertR[n];

    int edgeL[m];
    int edgeR[m];


    vertL[0]=1;
    vertR[0] = 0;
    struct edgeListNode* eln = (struct edgeListNode*) vertices[0]->outL;
    while(eln!=NULL){
        vertR[0] = vertR[0] + weightEdge[eln->ed->edgNum];
        eln = eln->next;
    }
    for(int i=1;i<n;i++){
        vertL[i] = vertR[0]; //maximum possible
        vertR[i] = 1;   //minimum possible
    }

    //Determining the leftmost and Right most chain passing through edges ad vertices
    for(int i=0;i<n;i++){
        int curr = vertL[i];
        struct edgeListNode* ol = vertices[i]->outL;
        while(ol!=NULL){
            int eno = ol->ed->edgNum;
            edgeL[eno] = curr;
            edgeR[eno] = curr + weightEdge[eno]-1;
            int p2num = ol->ed->p2->Ival;
            vertL[p2num] = vertL[p2num]<curr?vertL[p2num]:curr;
            vertR[p2num] = vertR[p2num]>(curr + weightEdge[eno]-1)?vertR[p2num]:(curr + weightEdge[eno]-1);
            curr = curr + weightEdge[eno];
            ol=ol->next;
        }
        //printf("Vertex Chain: %d %d\n",vertL[i],vertR[i]);
    }

    /*
    for(int i=0;i<m;i++){
        printf("Edge Chain: %d %d\n",edgeL[i],edgeR[i]);
    }
    */
    numchains = vertR[0];

    //Using SuperTree to store one edge w.r.t one chain
    struct superTreeNode* rootST = recConstrSuperTree(1,numchains,NULL);

    int edgesInChain[numchains+1];
    for(int i=0;i<=numchains;i++){
        edgesInChain[i]=0;
    }
    int edgeInSuperNode[m];
    for(int i=0;i<n;i++){
        struct edgeListNode* olh = vertices[i]->outL;
        while(olh!=NULL){
            int edn = olh->ed->edgNum;
            int leftv = edgeL[edn];
            int rightv = edgeR[edn];
            struct superTreeNode* lca = findLCA(rootST,leftv,rightv); //finding the chain of the edge that is LCA of all the chains passing through that edge
            int chn = lca->cindex;
            //printf("Chain %d\n",chn);
            edgesInChain[chn]++; //Number of edges stored in supertreenode for the given chain
            edgeInSuperNode[edn] = chn; //The chain in which the given edge is stored
            olh = olh->next;
        }
    }

    /*
    for(int i=0;i<=numchains;i++){
        printf("%d\t",edgesInChain[i]);
    }
    printf("\n");
    for(int i=0;i<m;i++){
        printf("%d\t",edgeInSuperNode[i]);
    }
    printf("\n");
    */
    initializeChainNodes(rootST,edgesInChain);
    int currVal[numchains+1]; //The current index to be used in the chain node for storing the edge
    for(int i=0;i<=numchains;i++){
        currVal[i] = 0;
    }
    for(int i=0;i<n;i++){
        struct edgeListNode* olh = vertices[i]->outL;
        while(olh!=NULL){
            struct edge* edgptr = olh->ed;
            int chainNum = edgeInSuperNode[edgptr->edgNum];
            assignToChainNum(rootST,edgptr,chainNum,currVal[chainNum]);
            currVal[chainNum]++;
            olh = olh->next;
        }
    }

    double xmin = vertices[0]->xcoord;
    double xmax = vertices[0]->xcoord;

    for(int i=1;i<n;i++){
        xmin = (vertices[i]->xcoord<xmin?vertices[i]->xcoord:xmin);
        xmax = (vertices[i]->xcoord>xmax?vertices[i]->xcoord:xmax);
    }
    //travSuperTree(rootST,edgesInChain);

    //query count
    int q;
    scanf("%d",&q);


    //Querying the points using super tree
    for(int qt=0;qt<q;qt++){
        double xq,yq;
        //input query point
        scanf("%lf",&xq);
        scanf("%lf",&yq);
        if(yq<vertices[0]->ycoord){
            printf("Query point is below the entire planar subdivision\n");
            continue;
        }
        else if(yq>vertices[n-1]->ycoord){
            printf("Query point is above the entire planar subdivision\n");
            continue;
        }
        else if(xq<xmin){
            printf("Query point is left to the entire planar subdivision\n");
            continue;
        }
        else if(xq>xmax){
            printf("Query point is right to the entire planar subdivision\n");
            continue;
        }

        struct superTreeNode* g = rootST;
        int l=0,r=numchains+1;
        int r1=0,r2=-1;
        //r1 -1 for left to 0 for on and +1 for right to segment
        //r2 edge number to which location relates
        while(1){
            int j = g->cindex;
            //printf("%d Subtree Node Entered",j);
            if(l>=j){
                g = g->right;
                continue;
            }
            if(r<=j){
                g=g->left;
                continue;
            }
            //printf("Processing chain %d\n",j);
            struct edge **gcNode = g->chNode;
            int numedges = edgesInChain[j];

            struct edge *projej=NULL;
            int le=0,re = numedges;
            while(le<=re){
                int mid = (le+re)/2;
                struct edge* ed = gcNode[mid];
                if(ed->p1->ycoord<=yq && ed->p2->ycoord>=yq){
                    projej = ed;
                    break;
                }
                else if(ed->p1->ycoord>yq){
                    re = mid-1;
                }
                else{
                    le = mid+1;
                }
            }


            double x1 = projej->p1->xcoord, y1 = projej->p1->ycoord;
            double x2 = projej->p2->xcoord, y2 = projej->p2->ycoord;
            //printf("For chain %d at %d: %lf %lf %lf %lf\n",j,projej->edgNum,x1,y1,x2,y2);

            if(y1==y2){
                if(x1>=xq && x2<=xq){
                    r1 = 0, r2 = projej->edgNum;
                    break;
                }

                if(x1<xq){
                    r1 = 1, r2 = projej->edgNum;
                    l = edgeR[r2];
                }
                else{
                    r1 = -1, r2 = projej->edgNum;
                    r = edgeL[r2];
                }

            }else{
                double d = (xq-x1)*(y2-y1)-(yq-y1)*(x2-x1);
                if(d<0){
                    r1=-1,r2=projej->edgNum;
                    r = edgeL[r2];
                }
                else if(d==0){
                    r1 = 0, r2 = projej->edgNum;
                    break;
                }
                else{
                    r1 = 1, r2 = projej->edgNum;
                    l = edgeR[r2];
                }
            }

            if(r-l<=1){
                break;
            }

        }

        double p1x = edges[r2].p1->xcoord,p1y = edges[r2].p1->ycoord;
        double p2x = edges[r2].p2->xcoord,p2y = edges[r2].p2->ycoord;
        //printf("Location: %d %d\n",r1,r2);
        printf("Query point is just %s the segment ((%lf,%lf),(%lf,%lf))\n",r1==1?"right of":(r1==0?"on":"left of"),p1x,p1y,p2x,p2y);

    }

}
