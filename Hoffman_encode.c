#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXVALUE  1000          //�����ı�����ַ�����
#define MAXLEAF   256           //���Ҷ�������������ͬ�ַ�����
#define MAXBIT    MAXLEAF-1     //������󳤶�
#define MAXNODE   MAXLEAF*2-1   //��������

typedef struct{        //Huffman�����ṹ��
    float weight;      //���Ȩֵ���������ַ����ֵ�Ƶ�ʣ���Ƶ��/�ַ�������
    int parent;        //�����λ����������ʼ-1
    int lchild;        //����λ����������ʼ-1
    int rchild;        //�Һ���λ����������ʼ-1
} HNodeType;

typedef struct{        //Huffman����ṹ��
    int bit[MAXBIT];   //�ַ��Ĺ���������
    int start;         //�ñ���������bit�еĿ�ʼλ��
} HCodeType;

int TextStatistics(char text[],float weight[]) {
    //ͳ��ÿ���ַ��ĳ���Ƶ�Σ����س��ֵĲ�ͬ�ַ��ĸ���
    //���ֵ��ַ������ch�У���Ӧ�ַ��ĳ���Ƶ�δ����weight��
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
	//����Ƶ������Ƶ��
	int index=0;
	while(weight[index]!=0){
		weight[index]/=strlen(text);
		index++;
	}
	//���� ch_index��ֵ��Ϊtext�ַ����в�ͬ�ַ��ĸ���
	return 5;
}

// �� HuffNodes[0..range]�У��ҵ���С�Ľ����������s1,s2 ���Ѿ��ҵ����Ľ��������������out[]��
void select(HNodeType HuffNodes[],int range,int *s1,int *s2){
	//���ҵ�һ����Сֵ ��
	float min1 = 5;
    int index1;
	for(index1=0;index1<=range;index1++){

		if(HuffNodes[index1].weight < min1 && HuffNodes[index1].parent ==-1){
			//�жϸý���Ƿ�ѡ��������ý��parentΪ0������Ϊ��ѡ
				min1 = HuffNodes[index1].weight;
				*s1 = index1 ;
		}
	}


	//�ҵ�2����Сֵ
	float min2 = 5;
	int index2;
	for(index2=0;index2<=range ;index2++){
		if(HuffNodes[index2].weight < min2 && HuffNodes[index2].parent ==-1 && index2!=*s1){
			//�жϸý���Ƿ�ѡ������Ҫ�ж����Ƿ�s1ѡ��
				min2 = HuffNodes[index2].weight;
				*s2 = index2 ;
		}
	}

}

//����һ��Huffman�������������HuffNodes��
void HuffmanTree(HNodeType HuffNodes[], float weight[], int n){

    if(n>MAXLEAF) {
    	printf("����Ҷ����������!\n");
    	return;
	}
	if(n<=1) return;

	int m = 2*n-1;//����ܸ���

	int node_index = 0;
	//�����Ҷ�ڵ�
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
	//�����Ҷ�ڵ�
	for(;node_index<m;node_index++){
		HuffNodes[node_index].weight = 0;
		HuffNodes[node_index].parent = -1;
		HuffNodes[node_index].lchild = -1;
		HuffNodes[node_index].rchild = -1;
	}

	//����Huffmantree

	int s1,s2,i;//��Сֵ����

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
    //����Huffman���룬Huffman��������HuffCodes��
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
    //Huffman�������������������������result�У�������һ��resultλ������
    //���ڵ� Ϊroot

	if (root!=-1){
		resultIndex = MidOrderTraverse( HuffNodes,result,HuffNodes[root].lchild,resultIndex);
		result[resultIndex++] = HuffNodes[root].weight;
		resultIndex = MidOrderTraverse( HuffNodes,result,HuffNodes[root].rchild,resultIndex);
	}

	return resultIndex;
}
*/
int main(){

    HNodeType HuffNodes[MAXNODE];   // ����һ�����ṹ������
    HCodeType HuffCodes[MAXLEAF];   // ����һ������ṹ������
    char text[MAXVALUE+1], ch[]="ATGCN";
    float weight[MAXLEAF], result[MAXNODE];
    int i, j, k, n, resultIndex;
    FILE *fp1,*fp2;
    fp1=fopen("randomSeq.txt","r");//��������ļ�
    fp2=fopen("huffmanSeq.txt","w+");

    while(fgets(text,1000,fp1)!=NULL)
    {
        if(text[0]!='>')
        {
            //�ַ�����n
            n = TextStatistics(text, weight);

            // �������������
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
    // ���Huffman��������������
    resultIndex = MidOrderTraverse(HuffNodes, result, 2*n-2, 0);
    printf("\nHuffman���������������ǣ�");

    for (i=0; i<resultIndex; i++)
        if (i < resultIndex-1)
            printf("%.4f, ", result[i]);
        else
            printf("%.4f\n\n", result[i]);
    */
    return 0;
}
