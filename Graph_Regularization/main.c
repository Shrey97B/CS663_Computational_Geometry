//Author Shrey Bhatt 20111060
#include <stdio.h>
#include <stdlib.h>
#include<math.h>

int N,M;

struct point{
    double xcoord;
    double ycoord;
    int iVal;
    struct edgeListNode* outL;
    struct edgeListNode* outT;
    struct edgeListNode* inL;
    struct edgeListNode* inT;
    int countOut;
    int countIn;
};

struct edgeListNode{
    struct edge* edg;
    struct edgeListNode* next;
};

struct edge{
    struct point* p1;
    struct point* p2;
    int eNum;
};

struct edgePointPairNode{
    struct edge* ed;
    struct point* p;
    int height;
    int count;
    struct edgePointPairNode* left;
    struct edgePointPairNode* right;
};

struct point** vertices;
struct edge* edges;

int pointcomp(const void* point1, const void* point2){
    const struct point **ver1 = (const struct point **) point1;
    const struct point **ver2 = (const struct point **) point2;
    double x1=(*ver1)->xcoord;
    double x2=(*ver2)->xcoord;
    double y1=(*ver1)->ycoord;
    double y2=(*ver2)->ycoord;

    if(y1!=y2){
        if(y1>y2){
            return 1;
        }
        return -1;
    }
    else{
        if(x2>x1){
            return 1;
        }
        return -1;
    }

}

void insertIntoOutgoing(struct point* p,struct edge* ed){
    struct edgeListNode* edH = p->outL;
    struct edgeListNode* edT = p->outT;

    if(edH==NULL){
        edH = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
        edH->edg = ed;
        edH->next = NULL;
        edT = edH;
    }
    else{
        struct edgeListNode* edn = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
        edn->edg = ed;
        edn->next = NULL;
        edT->next =edn;
        edT = edn;
    }

    p->outL = edH;
    p->outT = edT;
    p->countOut = p->countOut + 1;

}

void insertIntoIncoming(struct point* p,struct edge* ed){
    struct edgeListNode* edH = p->inL;
    struct edgeListNode* edT = p->inT;

    if(edH==NULL){
        edH = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
        edH->edg = ed;
        edH->next = NULL;
        edT = edH;
    }
    else{
        struct edgeListNode* edn = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
        edn->edg = ed;
        edn->next = NULL;
        edT->next =edn;
        edT = edn;
    }

    p->inL = edH;
    p->inT = edT;
    p->countIn = p->countIn + 1;

}

double getATan2Val(struct point* p1,struct point* p2){
    double ydiff = p2->ycoord-p1->ycoord;
    double xdiff = p2->xcoord-p1->xcoord;
    return atan2(ydiff,xdiff);
}

void spList(struct edgeListNode* head,
                    struct edgeListNode** lef, struct edgeListNode** rig)
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
    *lef = head;
    *rig = slow->next;
    slow->next = NULL;
}

struct edgeListNode* mergeOutGList(struct edgeListNode* lef, struct edgeListNode* rig)
{
    struct edgeListNode* res = NULL;

    if (lef == NULL)
        return (rig);
    else if (rig == NULL)
        return (lef);

    struct edge* lefed = lef->edg;
    struct edge* riged = rig->edg;

    double leftan2v = getATan2Val(lefed->p1,lefed->p2);
    double rigtan2v = getATan2Val(riged->p1,riged->p2);
    if (leftan2v >= rigtan2v) {
        res = lef;
        res->next = mergeOutGList(lef->next, rig);
    }
    else {
        res = rig;
        res->next = mergeOutGList(lef, rig->next);
    }
    return (res);
}

struct edgeListNode* mergeInGList(struct edgeListNode* lef, struct edgeListNode* rig)
{
    struct edgeListNode* res = NULL;

    if (lef == NULL)
        return (rig);
    else if (rig == NULL)
        return (lef);

    struct edge* lefed = lef->edg;
    struct edge* riged = rig->edg;

    double leftan2v = getATan2Val(lefed->p1,lefed->p2);
    double rigtan2v = getATan2Val(riged->p1,riged->p2);
    if (leftan2v <= rigtan2v) {
        res = lef;
        res->next = mergeInGList(lef->next, rig);
    }
    else {
        res = rig;
        res->next = mergeInGList(lef, rig->next);
    }
    return (res);
}

void sortOutGList(struct edgeListNode** refh){
    struct edgeListNode* head = *refh;
    struct edgeListNode* lef;
    struct edgeListNode* rig;

    if ((head == NULL) || (head->next == NULL)) {
        return;
    }

    spList(head, &lef, &rig);

    sortOutGList(&lef);
    sortOutGList(&rig);

    *refh = mergeOutGList(lef, rig);
}

void sortInGList(struct edgeListNode** refh){
    struct edgeListNode* head = *refh;
    struct edgeListNode* lef;
    struct edgeListNode* rig;

    if ((head == NULL) || (head->next == NULL)) {
        return;
    }

    spList(head, &lef, &rig);

    sortInGList(&lef);
    sortInGList(&rig);

    *refh = mergeInGList(lef, rig);
}

int rightOfOrOn(struct point* p, struct edgePointPairNode* epp){
    double xp = p->xcoord, yp = p->ycoord;
    double x1 = epp->ed->p1->xcoord, y1=epp->ed->p1->ycoord;
    double x2 = epp->ed->p2->xcoord, y2=epp->ed->p2->ycoord;
    double d = (xp-x1)*(y2-y1)-(yp-y1)*(x2-x1);

    if(d>=0){
        return 1;
    }
    else{
        return 0;
    }

}

struct edgePointPairNode* newEdgePointPairNode(struct edge* ed,struct point* p){
    struct edgePointPairNode* node = (struct edgePointPairNode*) malloc(sizeof(struct edgePointPairNode));
    node->count=1;
    node->height = 1;
    node->left = NULL;
    node->right = NULL;
    node->ed = ed;
    node->p=p;
    return node;
};

int height(struct edgePointPairNode* root){
    if(root==NULL){
        return 0;
    }
    return root->height;
}

int getBalance(struct edgePointPairNode *root)
{
    if (root == NULL)
        return 0;
    return height(root->left) - height(root->right);
}

int getCount(struct edgePointPairNode* root){
    if(root==NULL){
        return 0;
    }
    return root->count;
}

int maxn(int x1,int x2){
    if(x1>x2){
        return x1;
    }
    return x2;
}

struct edgePointPairNode *rightRotate(struct edgePointPairNode *y)
{
    struct edgePointPairNode *x = y->left;
    struct edgePointPairNode *T2 = x->right;

    x->right = y;
    y->left = T2;

    y->height = 1 + maxn(height(y->left),height(y->right));
    x->height = 1 + maxn(height(x->left),height(x->right));
    y->count = 1 + getCount(y->left) + getCount(y->right);
    x->count = 1 + getCount(x->left) + getCount(x->right);

    return x;
}

struct edgePointPairNode *leftRotate(struct edgePointPairNode *x)
{
    struct edgePointPairNode *y = x->right;
    struct edgePointPairNode *T2 = y->left;

    y->left = x;
    x->right = T2;

    x->height = 1 + maxn(height(x->left),height(x->right));
    y->height = 1 + maxn(height(y->left),height(y->right));
    x->count = 1 + getCount(x->left) + getCount(x->right);
    y->count = 1 + getCount(y->left) + getCount(y->right);

    return y;
}

struct edgePointPairNode* insertIntoAVLTree(struct edgePointPairNode* root, struct edge* ed, struct point* p,int nodeNum){

    if (root == NULL)
        return(newEdgePointPairNode(ed,p));

    int leftcount = getCount(root->left);

    if (nodeNum <= (leftcount+1))
        root->left  = insertIntoAVLTree(root->left, ed,p,nodeNum);
    else
        root->right = insertIntoAVLTree(root->right, ed,p,nodeNum - (leftcount+1));

    root->height = 1 + maxn(height(root->left),height(root->right));
    root->count = 1 + getCount(root->left) + getCount(root->right);

    int balance = getBalance(root);

    if (balance > 1 && nodeNum <= getCount(root->left->left))
        return rightRotate(root);
    else if (balance>1){
        root->left =  leftRotate(root->left);
        return rightRotate(root);
    }

    if (balance < -1 && (nodeNum - (getCount(root->left) + 1)) > (getCount(root->right->left) + 1))
        return leftRotate(root);
    else if (balance < -1 )
    {
        root->right = rightRotate(root->right);
        return leftRotate(root);
    }

    return root;
};

struct edgePointPairNode* minValueNode(struct edgePointPairNode* root)
{
    struct edgePointPairNode* curr = root;

    while (curr->left != NULL)
        curr = curr->left;

    return curr;
}

struct edgePointPairNode* deleteFromAVLTree(struct edgePointPairNode* root, int nodeNum, int nodeDelete)
{

    if (root == NULL)
        return root;

    if ( nodeDelete < nodeNum ){
        if(root->left!=NULL){
            root->left = deleteFromAVLTree(root->left, nodeNum - getCount(root->left->right) - 1, nodeDelete);
        }
    }
    else if( nodeDelete > nodeNum ){
        if(root->right!=NULL){
            root->right = deleteFromAVLTree(root->right, nodeNum + getCount(root->right->left) + 1, nodeDelete);
        }
    }
    else
    {
        if( (root->left == NULL) || (root->right == NULL) )
        {
            struct edgePointPairNode *temp = root->left ? root->left :
                                             root->right;

            if (temp == NULL)
            {
                temp = root;
                root = NULL;
            }
            else{
                *root = *temp;
            }
            free(temp);
        }
        else
        {
            struct edgePointPairNode* temp = minValueNode(root->right);
            root->ed = temp->ed;
            root->p = temp->p;
            root->right = deleteFromAVLTree(root->right, nodeNum + getCount(root->right->left) + 1,nodeDelete+1);
        }
    }

    if (root == NULL)
      return root;

    root->height = 1 + maxn(height(root->left),height(root->right));
    root->count = 1 + getCount(root->left) + getCount(root->right);

    int balance = getBalance(root);

    if (balance > 1 && getBalance(root->left) >= 0)
        return rightRotate(root);

    if (balance > 1 && getBalance(root->left) < 0)
    {
        root->left =  leftRotate(root->left);
        return rightRotate(root);
    }

    if (balance < -1 && getBalance(root->right) <= 0)
        return leftRotate(root);

    if (balance < -1 && getBalance(root->right) > 0)
    {
        root->right = rightRotate(root->right);
        return leftRotate(root);
    }

    return root;
}

struct point* getPforJVal(struct edgePointPairNode* root,int nodeNP,int nodeNum){
    if(root==NULL){
        return NULL;
    }

    if(nodeNP==nodeNum){
        return root->p;
    }

    if(nodeNP<nodeNum){
        if(root->left!=NULL){
            return getPforJVal(root->left,nodeNP,nodeNum - getCount(root->left->right) - 1);
        }
    }
    else{
        if(root->right!=NULL){
            return getPforJVal(root->right,nodeNP,nodeNum + getCount(root->right->left) + 1);
        }
    }

    return NULL;

}

void getJVal(struct edgePointPairNode* root,struct point* p,int* j,int nodeIndex){

    if(root==NULL){
        return;
    }

    if(root->ed->eNum==0 || rightOfOrOn(p,root)){
        *j = nodeIndex;
        if(root->right!=NULL){
            getJVal(root->right,p,j,nodeIndex + getCount(root->right->left)+1);
        }
    }
    else{
        if(root->left!=NULL){
            getJVal(root->left,p,j,nodeIndex - getCount(root->left->right) - 1);
        }
    }

}

void updateSequencePoint(struct edgePointPairNode* root,struct point* p,int nodeNum,int nodeUpdate){
    if(root==NULL){
        return;
    }

    if(nodeNum==nodeUpdate){
        root->p = p;
        return;
    }

    if(nodeUpdate<nodeNum){
        if(root->left!=NULL){
            updateSequencePoint(root->left,p,nodeNum - getCount(root->left->right) - 1,nodeUpdate);
        }
        return;
    }
    else{
        if(root->right!=NULL){
            updateSequencePoint(root->right,p,nodeNum + getCount(root->right->left) + 1,nodeUpdate);
        }
        return;
    }

}

void preOrderTraversal(struct edgePointPairNode* root){
    if(root==NULL){
        return;
    }

    //printf("%d %d %d %d\n",root->p->iVal,root->ed->p1->iVal,root->ed->p2->iVal,height(root));
    preOrderTraversal(root->left);
    preOrderTraversal(root->right);

}

int main(){

//input the number of vertices
    scanf("%d",&N);
    vertices = (struct point**) malloc((N+1)*sizeof(struct point*));
    for(int i=1;i<=N;i++){
        vertices[i] = (struct point *) malloc(sizeof(struct point));
        scanf("%lf",&(vertices[i]->xcoord));    //x coordinate of the vertex
        scanf("%lf",&(vertices[i]->ycoord));    //y coordinate of the vertex
        vertices[i]->outL = NULL;
        vertices[i]->inL = NULL;
        vertices[i]->outT = NULL;
        vertices[i]->inT = NULL;
        vertices[i]->countOut = 0;
        vertices[i]->countIn = 0;
    }

    scanf("%d",&M);
    edges = (struct edge*) malloc((M + 2*N) * sizeof(struct edge)); // (M+2N) because maximally 2N edges can be added after regularization
    int orig_m = M, pass1_add = 0,pass2_add = 0;

    for(int i=0;i<M;i++){
        int p1num,p2num;    //edges are represented by two 1-based indexes of above points
        scanf("%d",&p1num);
        scanf("%d",&p2num);
        edges[i].p1 = vertices[p1num];
        edges[i].p2 = vertices[p2num];
        edges[i].eNum = i+1;
    }

    qsort((void*)(vertices + 1),N,sizeof(vertices[1]),pointcomp);   //sorting the vertex according to the y-coordinates lower to high
    for(int i=1;i<=N;i++){
        vertices[i]->iVal = i;
    }

    for(int i=0;i<M;i++){
        int p1num = edges[i].p1->iVal;
        int p2num = edges[i].p2->iVal;
        if(p1num>p2num){    //keeping p1 as lower ordered vertex and p2 as higher one
            struct point* temp = vertices[p1num];
            edges[i].p1 = vertices[p2num];
            edges[i].p2 = temp;
        }
    }


    //construct incoming and outgoing edges
    for(int i=0;i<M;i++){
        insertIntoOutgoing(edges[i].p1,&edges[i]);
        insertIntoIncoming(edges[i].p2,&edges[i]);
    }

    for(int i=1;i<=N;i++){

        struct edgeListNode* oEdH = vertices[i]->outL;
        //sort incoming and outgoing edge list of vertices from left to right
        sortOutGList(&oEdH);
        struct edgeListNode* iEdH = vertices[i]->inL;
        sortInGList(&iEdH);
        vertices[i]->outL = oEdH;
        vertices[i]->inL = iEdH;
        if(oEdH!=NULL){
            do{
                vertices[i]->outT = oEdH;
                oEdH = oEdH->next;
            }while(oEdH!=NULL);
        }
        if(iEdH!=NULL){
            do{
                vertices[i]->inT = iEdH;
                iEdH = iEdH->next;
            }while(iEdH!=NULL);
        }
    }

    //trapezoids represented by pair of neighbouring edges enclosing it
    // x = - infinity edge to represent the open bounded trapezoid formed by it with the lef boundary of polygon
    struct edge* infEd = (struct edge*) malloc(sizeof(struct edge));
    struct point* infp1 = (struct point*) malloc(sizeof(struct point));
    struct point* infp2 = (struct point*) malloc(sizeof(struct point));
    infp1->xcoord = -INFINITY;
    infp1->ycoord = -INFINITY;
    infp1->iVal = -1;
    infp2->xcoord = -INFINITY;
    infp2->ycoord = INFINITY;
    infp2->iVal = -2;
    infEd->eNum=0;
    infEd->p1 = infp1;
    infEd->p2 = infp2;

    //Along with the edge also storing the point that will be at the top of trapezoid considering the area downwards 
    //This combination is stored in an AVL tree
    struct edgePointPairNode* trapSeqRoot = newEdgePointPairNode(infEd,vertices[N]);

    //First pass to add outgoing edges
    struct point* p0 = vertices[N];
    struct edgeListNode* iL = vertices[N]->inL;
    for(int i=1;i<=p0->countIn;i++){
        int lastcount = trapSeqRoot->count;
        trapSeqRoot = insertIntoAVLTree(trapSeqRoot,iL->edg,p0,lastcount+1);
        iL = iL->next;
    }

    //maintaining the count and list to store the new added edges
    int edgeAddCount = 0;
    struct edgeListNode* edgeAdded = NULL;
    struct edgeListNode* edgeATail = NULL;
    //printf("Traversal\n");
    //preOrderTraversal(trapSeqRoot);
    for(int i=N-1;i>=1;i--){
        p0 = vertices[i];
        int j=1;
        getJVal(trapSeqRoot,p0,&j,getCount(trapSeqRoot->left)+1);
        struct point* pconn = getPforJVal(trapSeqRoot,j,getCount(trapSeqRoot->left)+1);
        //printf("%d for point %d\n",j,i);

        //If there are no outgoing edges, the vertex is not regular and so need to add an outgoing edge
        if(p0->countOut==0){

            //printf("Add edge between %d and %d\n",p0->iVal,pconn->iVal);
            struct edge* ned = (struct edge*) malloc(sizeof(struct edge));
            ned->eNum = M + edgeAddCount;
            ned->p1 = vertices[i];
            ned->p2 = pconn;
            struct edgeListNode* nedn = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
            nedn->edg = ned;
            nedn->next = NULL;
            if(edgeAddCount==0){
                edgeAdded = edgeATail = nedn;
            }
            else{
                edgeATail->next = nedn;
                edgeATail = nedn;
            }
            edgeAddCount++;
        }

        updateSequencePoint(trapSeqRoot,p0,getCount(trapSeqRoot->left)+1,j-p0->countOut);

        struct edgeListNode* pass1_ol = p0->outL;
        for(int k=0;k<p0->countOut;k++){
            trapSeqRoot = deleteFromAVLTree(trapSeqRoot,getCount(trapSeqRoot->left) + 1,j-k);
        }

        j = j - p0->countOut;
        iL = p0->inL;
        for(int x = j+1;x<=j+p0->countIn;x++){
            //printf("Insert at %d\n",x);
            trapSeqRoot = insertIntoAVLTree(trapSeqRoot,iL->edg,p0,x);
            iL = iL->next;
        }

    }

    //Traversing the added edges into the edges array and incoming and outgoing list of vertices
    for(int i=M;i<M+edgeAddCount;i++){
        struct edge* curred = edgeAdded->edg;
        edges[i].eNum = i;
        edges[i].p1 = curred->p1;
        edges[i].p2 = curred->p2;
        edgeAdded = edgeAdded->next;
        insertIntoOutgoing(edges[i].p1,&edges[i]);
        insertIntoIncoming(edges[i].p2,&edges[i]);
        free(curred);

    }

    //sorting the incoming edges list after adding new vertices
    //no need to sort outgoing edges list
    for(int i=1;i<=N;i++){
        struct edgeListNode* iEdH = vertices[i]->inL;
        sortInGList(&iEdH);
        vertices[i]->inL = iEdH;
        if(iEdH!=NULL){
            do{
                vertices[i]->inT = iEdH;
                iEdH = iEdH->next;
            }while(iEdH!=NULL);
        }
    }

    M+=edgeAddCount;
    pass1_add=edgeAddCount;

    //initializing the tree for pass2 - incoming edges
    trapSeqRoot = newEdgePointPairNode(infEd,vertices[1]);

    p0 = vertices[1];
    struct edgeListNode* oL = vertices[1]->outL;
    for(int i=1;i<=p0->countOut;i++){
        int lastcount = trapSeqRoot->count;
        trapSeqRoot = insertIntoAVLTree(trapSeqRoot,oL->edg,p0,lastcount+1);
        oL = oL->next;
    }

    //initializing the below variables again for incoming edges
    edgeAddCount = 0;
    edgeAdded = NULL;
    edgeATail = NULL;
    //printf("Traversal\n");
    //preOrderTraversal(trapSeqRoot);
    for(int i=2;i<=N;i++){
        p0 = vertices[i];
        int j=1;
        getJVal(trapSeqRoot,p0,&j,getCount(trapSeqRoot->left)+1);
        struct point* pconn = getPforJVal(trapSeqRoot,j,getCount(trapSeqRoot->left)+1);
        //printf("%d for point %d\n",j,i);

        //Adding edges if number of incoming edges are zero
        if(p0->countIn==0){

            //printf("Add edge between %d and %d\n",p0->iVal,pconn->iVal);
            struct edge* ned = (struct edge*) malloc(sizeof(struct edge));
            ned->eNum = M + edgeAddCount;
            ned->p2 = vertices[i];
            ned->p1 = pconn;
            struct edgeListNode* nedn = (struct edgeListNode*) malloc(sizeof(struct edgeListNode));
            nedn->edg = ned;
            nedn->next = NULL;
            if(edgeAddCount==0){
                edgeAdded = edgeATail = nedn;
            }
            else{
                edgeATail->next = nedn;
                edgeATail = nedn;
            }
            edgeAddCount++;
        }

        updateSequencePoint(trapSeqRoot,p0,getCount(trapSeqRoot->left)+1,j-p0->countIn);

        struct edgeListNode* pass1_il = p0->inL;
        for(int k=0;k<p0->countIn;k++){
            trapSeqRoot = deleteFromAVLTree(trapSeqRoot,getCount(trapSeqRoot->left) + 1,j-k);
        }
        j = j - p0->countIn;
        oL = p0->outL;
        for(int x = j+1;x<=j+p0->countOut;x++){
            trapSeqRoot = insertIntoAVLTree(trapSeqRoot,oL->edg,p0,x);
            oL = oL->next;
        }

    }

    //Adding edges into edges array and incoming and outgoing edges of vertex for incoming edge pass
    for(int i=M;i<M+edgeAddCount;i++){
        struct edge* curred = edgeAdded->edg;
        edges[i].eNum = i;
        edges[i].p1 = curred->p1;
        edges[i].p2 = curred->p2;
        edgeAdded = edgeAdded->next;
        insertIntoOutgoing(edges[i].p1,&edges[i]);
        insertIntoIncoming(edges[i].p2,&edges[i]);
        free(curred);

    }

    for(int i=1;i<=N;i++){
        struct edgeListNode* oEdH = vertices[i]->outL;
        sortOutGList(&oEdH);
        vertices[i]->outL = oEdH;
        if(oEdH!=NULL){
            do{
                vertices[i]->outT = oEdH;
                oEdH = oEdH->next;
            }while(oEdH!=NULL);
        }
    }

    M+=edgeAddCount;
    pass1_add=edgeAddCount;

    if(M==orig_m){
        printf("No new edges have been added\n");
    }
    else{
        printf("New edges added between points:\n");
        for(int i=orig_m;i<M;i++){
            printf("(%lf,%lf) (%lf,%lf)\n",edges[i].p1->xcoord,edges[i].p1->ycoord,edges[i].p2->xcoord,edges[i].p2->ycoord);
        }
    }

    return 0;
}
