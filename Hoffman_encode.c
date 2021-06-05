#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXVALUE  1000          //输入文本最大字符个数
#define MAXLEAF   256           //最大叶结点个数，即最大不同字符个数
#define MAXBIT    MAXLEAF-1     //编码最大长度
#define MAXNODE   MAXLEAF*2-1   //最大结点个数

typedef struct{        //Huffman树结点结构体
    float weight;      //结点权值，这里是字符出现的频率，及频次/字符种类数
    int parent;        //父结点位置索引，初始-1
    int lchild;        //左孩子位置索引，初始-1
    int rchild;        //右孩子位置索引，初始-1
} HNodeType;

typedef struct{        //Huffman编码结构体
    int bit[MAXBIT];   //字符的哈夫曼编码
    int start;         //该编码在数组bit中的开始位置
} HCodeType;

int TextStatistics(char text[],float weight[]) {
    //统计每种字符的出现频次，返回出现的不同字符的个数
    //出现的字符存放在ch中，对应字符的出现频次存放在weight中
	int i;

	for(i=0;i<strlen(text);i++)
    {
        if(text[i]=='A'||text[i]=='a')
            weight[0]+=1;
        else if(text[i]=='T'||text[i]=='t')
            weight[1]+=1;
        else if(text[i]=='G'||text[i]=='g')
            weight[2]+=1;
        else if(text[i]=='C'||text[i]=='c')
            weight[3]+=1;
        else
            weight[4]+=1;
    }
	//根据频数计算频率
	int index=0;
	while(weight[index]!=0){
		weight[index]/=strlen(text);
		index++;
	}
	//最终 ch_index的值即为text字符串中不同字符的个数
	return 5;
}

// 从 HuffNodes[0..range]中，找到最小的结点索引赋给s1,s2 。已经找到过的结点索引被储存在out[]中
void select(HNodeType HuffNodes[],int range,int *s1,int *s2){
	//先找第一个最小值 。
	float min1 = 5;
    int index1;
	for(index1=0;index1<=range;index1++){

		if(HuffNodes[index1].weight < min1 && HuffNodes[index1].parent ==-1){
			//判断该结点是否被选过。如果该结点parent为0，则其为被选
				min1 = HuffNodes[index1].weight;
				*s1 = index1 ;
		}
	}


	//找第2个最小值
	float min2 = 5;
	int index2;
	for(index2=0;index2<=range ;index2++){
		if(HuffNodes[index2].weight < min2 && HuffNodes[index2].parent ==-1 && index2!=*s1){
			//判断该结点是否被选过。还要判断其是否被s1选了
				min2 = HuffNodes[index2].weight;
				*s2 = index2 ;
		}
	}

}

//构造一棵Huffman树，树结点存放在HuffNodes中
void HuffmanTree(HNodeType HuffNodes[], float weight[], int n){

    if(n>MAXLEAF) {
    	printf("超出叶结点最大数量!\n");
    	return;
	}
	if(n<=1) return;

	int m = 2*n-1;//结点总个数

	int node_index = 0;
	//构造各叶节点
	for(;node_index < n;node_index++){
	/*
		HuffNodes[node_index]->weight
		HuffNodes[node_index]->parent
		HuffNodes[node_index]->lchild
		HuffNodes[node_index]->rchild
	*/
		HuffNodes[node_index].weight = weight[node_index];
		HuffNodes[node_index].parent = -1;
		HuffNodes[node_index].lchild = -1;
		HuffNodes[node_index].rchild = -1;
	}
	//构造非叶节点
	for(;node_index<m;node_index++){
		HuffNodes[node_index].weight = 0;
		HuffNodes[node_index].parent = -1;
		HuffNodes[node_index].lchild = -1;
		HuffNodes[node_index].rchild = -1;
	}

	//构建Huffmantree

	int s1,s2,i;//最小值索引

	for(i = n;i < m;i++) {
		select(HuffNodes,i-1,&s1,&s2);
		HuffNodes[s1].parent = i;
		HuffNodes[s2].parent = i;
		HuffNodes[i].lchild = s1;
		HuffNodes[i].rchild = s2;
		HuffNodes[i].weight = HuffNodes[s1].weight + HuffNodes[s2].weight;

	}
}

void HuffmanCode(HNodeType HuffNodes[], HCodeType HuffCodes[], int n) {
    //生成Huffman编码，Huffman编码存放在HuffCodes中
    int start,i,c,f;
	for(i =0 ;i<n;i++){
		start = n-2;
		for(c = i , f=HuffNodes[i].parent ; f!=-1; c =f,f=HuffNodes[f].parent){
			if(c == HuffNodes[f].lchild) HuffCodes[i].bit[start--]=0;
			else HuffCodes[i].bit[start--]=1;
		}
		HuffCodes[i].start = start+1;
	}
}
/*
int MidOrderTraverse(HNodeType HuffNodes[], float result[], int root, int resultIndex) {
    //Huffman树的中序遍历，遍历结果存放在result中，返回下一个result位置索引
    //根节点 为root

	if (root!=-1){
		resultIndex = MidOrderTraverse( HuffNodes,result,HuffNodes[root].lchild,resultIndex);
		result[resultIndex++] = HuffNodes[root].weight;
		resultIndex = MidOrderTraverse( HuffNodes,result,HuffNodes[root].rchild,resultIndex);
	}

	return resultIndex;
}
*/
int main(){

    HNodeType HuffNodes[MAXNODE];   // 定义一个结点结构体数组
    HCodeType HuffCodes[MAXLEAF];   // 定义一个编码结构体数组
    char text[MAXVALUE+1], ch[]="ATGCN";
    float weight[MAXLEAF], result[MAXNODE];
    int i, j, k, n, resultIndex;
    FILE *fp1,*fp2;
    fp1=fopen("randomSeq.txt","r");//输入测试文件
    fp2=fopen("huffmanSeq.txt","w+");

    while(fgets(text,1000,fp1)!=NULL)
    {
        if(text[0]!='>')
        {
            //字符总数n
            n = TextStatistics(text, weight);

            // 输出哈夫曼编码
            HuffmanTree(HuffNodes, weight, n);
            HuffmanCode(HuffNodes, HuffCodes, n);

            for (i=0; i<n; i++) {
                fprintf(fp2,"%c", ch[i]);
                for(j=HuffCodes[i].start; j<n-1; j++)
                    fprintf(fp2,"%d", HuffCodes[i].bit[j]);
            }
            fprintf(fp2,"\n");
            for(i=0;i<strlen(text);i++)
            {
                if(text[i]=='a'||text[i]=='A')
                    for(j=HuffCodes[0].start; j<n-1; j++)
                        fprintf(fp2,"%d", HuffCodes[0].bit[j]);
                else if(text[i]=='t'||text[i]=='T')
                    for(j=HuffCodes[1].start; j<n-1; j++)
                        fprintf(fp2,"%d", HuffCodes[1].bit[j]);
                else if(text[i]=='g'||text[i]=='G')
                    for(j=HuffCodes[2].start; j<n-1; j++)
                        fprintf(fp2,"%d", HuffCodes[2].bit[j]);
                else if(text[i]=='c'||text[i]=='C')
                    for(j=HuffCodes[3].start; j<n-1; j++)
                        fprintf(fp2,"%d", HuffCodes[3].bit[j]);
                else
                    for(j=HuffCodes[4].start; j<n-1; j++)
                        fprintf(fp2,"%d", HuffCodes[4].bit[j]);
            }
            fprintf(fp2,"\n");
            for(k=0;k<n;k++)
                weight[k]=0;
        }
    }
    /*
    // 输出Huffman树的中序遍历结果
    resultIndex = MidOrderTraverse(HuffNodes, result, 2*n-2, 0);
    printf("\nHuffman树的中序遍历结果是：");

    for (i=0; i<resultIndex; i++)
        if (i < resultIndex-1)
            printf("%.4f, ", result[i]);
        else
            printf("%.4f\n\n", result[i]);
    */
    return 0;
}
