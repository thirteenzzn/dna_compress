#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int n=4;  //默认匹配长度

typedef struct PolyNode *Polynomial;
//存储DNA序列，用换行符来区分前后不同序列d
struct PolyNode{
    Polynomial pre;  //左结点
    char base;  //碱基
    int num;  //位置
    Polynomial next;  //右结点
};









//10进制整型数转换成字符串d
void Int2String(int num,char *str,int *digit)
{
	int i = 0;
	if(num<0)//如果num为负数，将num变正 
	{
		num = -num;
		str[i++] = '-';
	} 
	//转换 
	do
	{
		str[i++] = num%10+'0';
		num /= 10;	
	}while(num);
	
	str[i] = '\0';
    *digit=i;
	
	//确定开始调整的位置 
	int j = 0;
	if(str[0]=='-')//如果有负号，负号不用调整 
	{
		j = 1;
		++i;
	}
	//对称交换 
	for(;j<i/2;j++)
	{ 
		str[j] = str[j] + str[i-1-j];
		str[i-1-j] = str[j] - str[i-1-j];
		str[j] = str[j] - str[i-1-j];
	} 
}









//把值加入链表中d
void Attach(char b, Polynomial *pp)
{
    Polynomial p;
    p=(Polynomial)malloc(sizeof(struct PolyNode));

    //存储序列
    if((*pp)->num==-1)
    {
        p->base=b;
        p->num=0;
        p->pre=*pp;
        p->next=NULL;
    }
    else
    {
        p->base=b;
        p->num=(*pp)->num+1;
        p->pre=*pp;
        p->next=NULL;
    }

    //修改pp的值
    (*pp)->next=p;
    *pp=p;
}









//读取DNA序列d
Polynomial ReadFile(char *filename)
{
    Polynomial p,prear,t;
    FILE *fp;  //打开文件的指针
    FILE *fp1;  //存储DNA信息头的文件的指针
    FILE *fp2;  //存储DNA序列的文件的指针
    char str[80];
    int flag=0;

    p=(Polynomial)malloc(sizeof(struct PolyNode));
    p->next=NULL;
    p->pre=NULL;
    prear=p;
    fp=fopen(filename,"r");
    fp1=fopen("/home/aubin/bi296/project/ids.txt","w+");  //存储信息头
    fp2=fopen("/home/aubin/bi296/project/seq.txt","w+");  //存储序列
    while(fgets(str,80,fp)!=NULL)
    {
        if(str[0]=='>')  //判断是否为信息头
        {
            fwrite(str,sizeof(char),strlen(str)*sizeof(char),fp1);
            flag=1;
        }
        else
        {   
            if(str[strlen(str)-1]=='\n')
                str[strlen(str)-1]='\0';  //删除\n
            //存入链表
            int i;
            for(i=0;i<strlen(str);i++)
            {
                if(flag==1)
                {
                    flag=0;
                    prear->num=-1;  //初始化位置
                    Attach('\n',&prear);
                    Attach(str[i],&prear);
                }
                else
                {
                    Attach(str[i],&prear);  //将碱基插入链表尾部
                }
            }

            //写入文件
            fwrite(str,sizeof(char),strlen(str)*sizeof(char),fp2);
        }
    }
    //删除临时生成的头结点
    t=p;
    p=p->next;
    free(t);

    fclose(fp);
    fclose(fp1);
    fclose(fp2);
    return p;
}









//二进制压缩d
void b_compress(Polynomial p,int l)
{
    Polynomial ptmp;  //临时结点
    int i;
    char str[l*2];
    ptmp=p;
    for(i=0;i<l-1;i++)
    {
        if(i==0)
        {
            switch(ptmp->base)
            {
                case 'A': strcpy(str,"00"); break;
                case 'T': strcpy(str,"01"); break;
                case 'C': strcpy(str,"10"); break;
                case 'G': strcpy(str,"11"); break;
                case '\n': i=i;
            }
        }
        else
        {
            switch(ptmp->base)
            {
                case 'A': strcat(str,"00"); break;
                case 'T': strcat(str,"01"); break;
                case 'C': strcat(str,"10"); break;
                case 'G': strcat(str,"11"); break;
                case '\n': i=i;
            }
        }
        ptmp=ptmp->next;
    }

    //写入压缩文件
    FILE *fp;
    fp=fopen("/home/aubin/bi296/project/compress.txt","w+");
    fputs(str,fp);
    fclose(fp);
    free(ptmp);
}









//定位部分序列结尾d
void Findtail(Polynomial p,Polynomial *ptail)
{
    int i;
    *ptail=p;
    for(i=0;i<n-1;i++)
    {
        *ptail=(*ptail)->next;
    }
}









//寻找种子及子串d
int FindSubstr(Polynomial *pseed, Polynomial *pstr)
{
    int i;
    int flag=1;
    int key=1;  //判断是否找到子串
    Polynomial pseed1,pstr1,ptailseed,ptailstr,p;
    pstr1=*pstr;   
    ptailseed=*pseed;
    ptailstr=NULL;
    p=*pseed;
    Findtail(p,&ptailseed);

    while(ptailseed->next!=NULL)
    {   
        pseed1=p;
        Findtail(p,&pstr1);
        pstr1=pstr1->next;
        Findtail(pstr1,&ptailstr);
        while(ptailstr->next!=NULL)
        {
            flag=1;  //匹配正确
            for(i=0;i<n;i++)
            {
                if(pseed1->base!=pstr1->base)
                {
                    flag=0;  //匹配失败
                    break;
                }
                pseed1=pseed1->next;
                pstr1=pstr1->next;
            }
            if(flag==1)
                break;
            else
            {
                pstr1=pstr1->next;
                pseed1=p;
            }
            Findtail(pstr1,&ptailstr);
        }
        if(flag==1)
            break;
        if(ptailstr==NULL)
        {
            p=p->next;
            Findtail(p,&ptailseed);
        }
    }
    
    for(i=0;i<n;i++)
        pstr1=pstr1->pre;
    if(ptailseed->next==NULL)
    {
        key=0;
        return key;
    }

    *pstr=pstr1;

    //free(p);
    //free(pseed1);
    //free(pstr1);
    //free(ptailseed);
    //free(ptailstr);

    return key;
}









//删除子串d
void DeleteSeed(Polynomial *pstr,int l)
{
    Polynomial tmp;
    int i;    

    //删除子串
    tmp=*pstr;
    for(i=0;i<l;i++)
        tmp=tmp->next;
    (*pstr)->pre->next=tmp;
    tmp->pre=(*pstr)->pre;
    *pstr=tmp;
}









//比较子串与种子的不匹配个数
int compare(Polynomial pseedl,Polynomial pseedr,Polynomial pstrl,Polynomial pstrr)
{
    int length,i=0,sum=0;
    Polynomial seed,str;
    seed=pseedl;
    str=pstrl;
    length=pseedr->num-pseedl->num+1;
    for(i=0;i<length;i++)
    {
        if(seed->base!=str->base)
            sum++;
        seed=seed->next;
        str=str->next;
    }
    sum+=pstrr->num-pstrl->num+1-length;
    return sum;
}









//十进制转二进制
void TentoTwo(int x,char *s)
{
    int k=0;

    while(x)
    {
        if(x%2==0)
            *(s+k)='0';
        else
            *(s+k)='1';
        x=x/2;
        k++;
    }
    
    for(;k<4;k++)
    {
        *(s+k)='0';
    }
    *(s+k)='\0';
}









//比较字符串
void CompareStr(char *str1,char *str2,char *tumple1, char *tumple2)
{
    char *c;
    int i;
    c=strstr(str2,str1);

    if(c==str2)  //类型为插入，位置在后两位
    {
        char *s;
        s=(char *)malloc(sizeof(char)*5);
        TentoTwo(strlen(str2)-1,s);
        strcpy(tumple1,s);
        strcat(tumple1,"00");
        free(s);
        char *ss;
        ss=(char *)malloc(sizeof(char)*5);
        TentoTwo(strlen(str2),ss);
        strcpy(tumple2,ss);
        strcat(tumple2,"00");
        free(s);
    }
    else if(c==NULL)
    {
        for(i=0;i<strlen(str1);i++)
        {
            if(str1[i]!=str2[i])
                break;
        }
        char *s;
        s=(char *)malloc(sizeof(char)*5);
        TentoTwo(i,s);
        strcpy(tumple1,s);
        
        strcat(tumple1,"10");
        free(s);
        char *ss;
        ss=(char *)malloc(sizeof(char)*5);
        TentoTwo(strlen(str2),ss);
        strcpy(tumple2,ss);
        strcat(tumple2,"00");
        free(s);
    }
    else
    {
        char *s;
        s=(char *)malloc(sizeof(char)*5);
        TentoTwo(0,s);
        strcpy(tumple1,s);
        strcat(tumple1,"00");
        free(s);
        char *ss;
        ss=(char *)malloc(sizeof(char)*5);
        TentoTwo(1,ss);
        strcpy(tumple2,ss);
        strcat(tumple2,"00");
        free(s);
        
    }
}









//延伸并存入离线文件
void find(Polynomial p,Polynomial pseed, Polynomial *pstr)
{
    Polynomial pseedl,pstrl;  //向左延申的临时结点
    Polynomial pseedr,pstrr;  //向右延申的临时结点
    Polynomial pseedtail,pstrtail;  //种子、子串的结尾
    int i;
    int m=0;  //不匹配次数
    int flag=1;
    char tumple1[7]={'e','e','e','e','e','e','e'};  //信息元组1，初始化为e
    char tumple2[7]={'e','e','e','e','e','e','e'};  //信息元组2，初始化为e

    //向左延伸
    pseedl=pseed->pre;
    pstrl=(*pstr)->pre;
    //定位种子结尾
    Findtail(pseed,&pseedtail);
    //定位子串结尾
    Findtail(*pstr,&pstrtail);
    pseedr=pseedtail;
    pstrr=pstrtail;

    while(pseedl!=NULL&&pseedtail!=pstrl)  //排除末尾和种子与子串重叠的情况
    {
        if(pseedl->base=='\n')
        {
            pseedl=pseedl->next;
            pstrl=pstrl->next;
            break;
        }
        m=compare(pseedl,pseedr,pstrl,pstrr);
        if(m>=2)  //向左匹配结束
            break;

        //继续延伸
        pstrl=pstrl->pre;
        pseedl=pseedl->pre;
    }

    //向右延伸
    pseedr=pseedtail->next;
    pstrr=pstrtail->next;
    if(pstrr==NULL)
        pstrr=pstrr->pre;
    else if(pstrr->next==NULL)
    {
        if(m==0||m==1)
            pseedtail=pseedtail->next;
    }
    else
    {
        if(m==0)
            pseedtail=pseedtail->next->next;
        if(m==1)
            pseedtail=pseedtail->next;
    }
    
    while(pseedtail!=pstrl&&pstrr!=NULL)  //排除末尾和种子与子串重叠的情况
    {
        if(pstrr->base!=pseedr->base)  //不匹配
            m++;
        
        if(m>=2) //m大于2，终止
        {
            pseedr=pseedr->pre;
            pstrr=pstrr->pre;
            break;
        }
        pstrr=pstrr->next;
        pseedr=pseedr->next;
    }



    int l;  //子串长度
    l=pstrr->num-pstrl->num+1;
    
    //字典更新
    FILE *fp;
    fp=fopen("/home/aubin/bi296/project/dictionary.txt","a+");
    char *str;
    str=(char *)malloc(sizeof(char)*100);

    //存储种子
    Polynomial tmp;
    tmp=pseed;
    for(i=0;i<n;i++)
    {
        str[i]=tmp->base;
        tmp=tmp->next;
    }
    str[i]='\0';
    free(tmp);

    //格式换位及匹配类型(模糊匹配)
    strcat(str,"\tApprox\t");
    
    //存储位置
    char *tmpstr1;
    int digitn;
    tmpstr1=(char *)malloc(sizeof(char)*5);
    Int2String(pstrl->num,tmpstr1,&digitn);  //数字转换为字符串
    
    strcat(str,tmpstr1);  //存储位置
    strcat(str,"\t");  //格式控制


    //存储长度
    char *tmpstr;
    int digitl;
    tmpstr=(char *)malloc(sizeof(char)*5);
    Int2String(pstrr->num-pstrl->num+1,tmpstr,&digitl);  //数字转换为字符串

    strcat(str,tmpstr);  //存储长度
    strcat(str,"\t");  //格式控制

    //比对种子与子串

    /*换成字符串方便比较*/
    char *str1;
    char *str2;
    str1=(char *)malloc(sizeof(char)*(pseedtail->num-pseed->num+1));
    str2=(char *)malloc(sizeof(char)*(pstrr->num-pstrl->num+1));
    Polynomial t;
    t=pseedl;
    for(i=0;i<pseedtail->num-pseedl->num+1;i++)
    {
        str1[i]=t->base;
        t=t->next;
    }
    str1[i]='\0';
    free(t);

    Polynomial tt;
    tt=pstrl;
    for(i=0;i<pstrr->num-pstrl->num+1;i++)
    {
        str2[i]=tt->base;
        tt=tt->next;
    }
    str2[i]='\0';
    free(tt);
    /*比较字符串，存储信息元组*/
    CompareStr(str1,str2,tumple1,tumple2);

    //存储信息元组
    if(tumple1[0]=='e')
    {
        strcat(str,"\n"); 
    }
    else if(tumple2[0]=='e')
    {
        strcat(str,tumple1); 
        strcat(str,"\n"); 
    }
    else
    {
        strcat(str,tumple1); 
        strcat(str,"\t"); 
        strcat(str,tumple2); 
        strcat(str,"\n"); 
    }
    
    //写进文件
    fputs(str,fp);
    //关闭文件
    fclose(fp);
    //删除子串
    DeleteSeed(pstr,l);

    free(pseedl);
    free(pstrl);
    free(pseedr);
    free(pstrr);
    free(pseedtail);
    free(pstrtail);
}










//获取文件大小d
int file_size(char* filename)
{
    FILE *fp=fopen(filename,"r");
    if(!fp) return -1;
    fseek(fp,0L,SEEK_END);
    int size=ftell(fp);
    fclose(fp);
    
    return size;
}









//计算文件压缩比d
float fcomprese()
{
    int newnum=0,oldnum=0;  //初始化
    float rate=1.0;  //初始化比例
    newnum=file_size("/home/aubin/bi296/project/dictionary.txt")+file_size("/home/aubin/bi296/project/compress.txt")+file_size("/home/aubin/bi296/project/ids.txt");
    oldnum=file_size("/home/aubin/bi296/project/dna.fasta");
    rate=1.0*newnum/oldnum;
    return rate;
}









int main()
{
    Polynomial pd;  //存储序列
    Polynomial pseed;  //指向种子
    Polynomial pstr;  //指向子串
    Polynomial ptail;  //指向链表尾部
    Polynomial tmp;  //临时结点
    float srate,frate;  //压缩结果
    int length=0;  //剩余序列长度
    int key=1;
    int i=0;
    //读入文件
    char filename[80]="/home/aubin/bi296/project/dna.fasta";
    //scanf("%s",filename);

    pd=ReadFile(filename);  //读入序列
    //初始化种子和子串
    pseed=pd->next;
    Findtail(pstr,&ptail);
    while(ptail!=NULL)
    {
        pstr=pseed->next->next->next->next;
        key=FindSubstr(&pseed,&pstr);  //寻找种子及子串
        find(pd,pseed,&pstr);  //延伸并存入离线字典文件
        //printf("%d\n",key);
        Findtail(pstr,&ptail);
    }
        
    //计算剩余序列长度
    tmp=pd->next;
    while(tmp!=NULL)
    {
        length++;
        tmp=tmp->next;
    }
    free(tmp);
    
    b_compress(pd->next,length);  //二进制压缩剩余序列
    frate=fcompresse();  //计算文件压缩大小

    printf("%f\n",frate);
    return 0;
}