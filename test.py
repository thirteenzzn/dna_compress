import cv2
import streamlit as st
import numpy as np
import pandas as pd
from PIL import Image

image_bk=Image.open('5.jpg')
pse1=Image.open('1.jpg')
pse2=Image.open('2.png')
pse3=Image.open('3.png')
video_file=open('1.mp4','rb')
video=video_file.read()

st.image(image_bk,use_column_width=True)

"""
# DNA algorithm
"""

st.write('### Group3\n')

contentf='''

        展望
        '''

bk='''
    ## 背景
    DNA 测序在许多生物信息学应用中起着至关重要的作用，自 Sanger 测序方法被采用以来产生了重大影响。基于阵列的焦磷酸测序、合成测序和连接测序技术成倍地增加了 DNA 序列数据的生成。测序技术正在以比旧技术高得多的速度开发和产生测序数据。这种进步导致了越来越多的大型数据项目的经济发展。这些项目包括 1000 个基因组项目，国际癌症基因组项目，ENCODE 项目等。
    \n
    对 DNA 序列的分析也能对这些序列进行智能分析。压缩在高效的序列分类中也起着重要的作用。DNA 序列由 4 个碱基组成，A(腺嘌呤)，C(胞嘧啶)，G(鸟嘌呤)和 T(胸腺嘧啶)，两个比特足以代表每个碱基。此外，在 DNA 序列中发现的重复序列并不总是精确的；它们可以是不同的类型，如近似、反向、互补、反向互补和串联。而且这些重复的时间很长，频率也较低。传统的文本压缩算法只能有效地捕获短而频繁的重复；因此，用它们来压缩 DNA 序列往往会导致相同序列的扩增。所以，找到 DNA 序列中所有不同类型的重复序列并对它们进行编码以获得良好的压缩比是一项具有挑战性的任务
    \n
    ## 目的
    下一代基因测序技术迅猛发展为生命科学领域带来新生与活力，然而困境与挑战也随之而来。基因测序成本大幅降低，从而导致测序数据量迅速攀升，其增速已经远超计算机硬件成本的下降速度。基因测序产生的各类文件在存储和传输时由于数据量较大而消耗大量时间，从而导致出现明显的延迟，为后续数据的研究及使用带来诸多不便。可以预见，测序数据量的持续增长所带来的数据存储负担，将成为诸多研究机构，医疗机构甚至是个人无法回避的难题，因此研究针对基因测序数据的压缩算法是十分有必要的。
    '''

content_include1='''
#include <algorithm>
#include <condition_variable>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <utility>
#include <vector>

typedef unsigned char u8;
typedef unsigned short u16;
typedef std::map<u16, std::vector<int>> SegMap;
typedef std::pair<u16, std::vector<int>> SegPair;
typedef std::unique_lock<std::mutex> MutLock;

const int kThreadN = 8;
const int kMaxBuf = 8000;
const int kSegLen = 6;
const int kFindN = 20;
int threadN;
int stage;
std::mutex stage_mu;
std::condition_variable stage_cv;
int finishN;
std::mutex finishN_mu;
std::condition_variable finishN_cv;

char bufs[kThreadN][kMaxBuf + 1];
int bufs_len[kThreadN];
std::vector<SegPair> freqs[kThreadN];
u16 freq_result;

        '''
content_encode_one='''
static inline u16 encode_one(char ch) {
  switch (ch) {
    case 'A':
      return 0;
    case 'T':
      return 1;
    case 'C':
      return 2;
    case 'G':
      return 3;
    default:
      std::cerr << "[error]: find undefined char: " << ch << "\\n";
      return 0;
  }
}
'''
content_encode_foeward='''
static inline u16 encode_forward(int key, char ch) {
  key <<= 2;
  key += encode_one(ch);
  return key;
}
'''

content_encode='''
static u16 encode(char* seq) {
  u16 key = 0;
  for (int i = 0; i < kSegLen; i++) {
    key <<= 2;
    key += encode_one(seq[i]);
  }
  return key;
}
'''

content_adds_map='''
static inline void adds_map(SegMap& segmap, SegPair& result) {
  auto key = result.first;
  auto& locs = result.second;
  auto iter = segmap.find(key);
  if (iter == segmap.end()) {
    segmap[key] = std::move(locs);
  } else {
    auto& old_locs = iter->second;
    int old_last_loc = *old_locs.rbegin();
    int first_loc = locs[0];
    auto insert_begin = locs.begin();
    if (first_loc < old_last_loc) {
      insert_begin++;
    }
    old_locs.insert(old_locs.end(), insert_begin, locs.end());
  }
}
'''

content_add_map='''
static inline void add_map(SegMap& segmap, u16 key, int loc) {
  auto iter = segmap.find(key);
  if (iter == segmap.end()) {
    segmap[key] = { loc};
  } else {
    int last_loc = *iter->second.rbegin();
    // 避免两个segment重叠在一起
    if (loc >= last_loc + kSegLen) {
      iter->second.push_back(loc);
    }
  }
}
'''

content_struct='''
struct SegPairCmp {
  bool operator()(const SegPair& a, const SegPair& b) {
    return a.second.size() < b.second.size();
  }
};
'''

content_find_most_freq='''
static std::vector<SegPair> find_most_freq(SegMap& segmap, char* seq, int len) {
  u16 key = encode(seq);
  add_map(segmap, key, 0);
  for (int i = 1; i <= len - kSegLen; i++) {
    key = encode_forward(key, seq[i + kSegLen - 1]);
    add_map(segmap, key, i);
  }
  std::priority_queue<SegPair, std::vector<SegPair>, SegPairCmp> pq(
      std::make_move_iterator(segmap.begin()),
      std::make_move_iterator(segmap.end()));
  segmap.clear();
  auto max_ele = pq.top();
  pq.pop();
  int bound = std::max(int(max_ele.second.size() / kThreadN), 1);
  std::vector<SegPair> results;
  segmap.insert(max_ele);
  results.push_back(std::move(max_ele));
  while (!pq.empty() && pq.top().second.size() > bound) {
    segmap.insert(pq.top());
    results.push_back(std::move(pq.top()));
    pq.pop();
  }
  return results;
}
'''

content_remove_freq='''
static void remove_freq(const std::vector<int>& locs, char* seq, int& len) {
  int prev_loc = -kSegLen;
  char *p = seq, *q = seq;
  for (int loc : locs) {
    int cp_len = loc - prev_loc - kSegLen;
    memcpy(q, p, cp_len);
    q += cp_len;
    p += cp_len + kSegLen;
    prev_loc = loc;
  }
  memcpy(q, p, len - prev_loc - kSegLen);
  len -= locs.size() * kSegLen;
}
'''
content_thread_routine='''
static void thread_routine(int tid) {
  auto& buf = bufs[tid];
  auto& buf_len = bufs_len[tid];
  auto& freq = freqs[tid];
  bool quit = false;
  MutLock stage_lck(stage_mu, std::defer_lock);
  MutLock finishN_lck(finishN_mu, std::defer_lock);
  SegMap segmap;
  int cur_stage = 1;
  while (true) {
    finishN_lck.lock();
    finishN++;
    finishN_cv.notify_one();
    finishN_lck.unlock();
    if (quit) break;
    stage_lck.lock();
    stage_cv.wait(stage_lck, [&quit, cur_stage]() {
      if (stage == -1) quit = true;
      return stage == -1 || stage == cur_stage;
    });
    cur_stage++;
    stage_lck.unlock();
    if (quit) {
      int aaa = 1;
    }
    auto iter = segmap.find(freq_result);
    if (iter != segmap.end()) {
      remove_freq(iter->second, buf, buf_len);
    }
    if (!quit) {
      segmap.clear();
      auto comp_results = find_most_freq(segmap, buf, buf_len);
      int offset = 0;
      for (int i = 0; i < tid; i++) {
        offset += bufs_len[i];
      }
      for (auto& result : comp_results) {
        for (auto& loc : result.second) {
          loc += offset;
        }
      }
      freq = std::move(comp_results);
    }
  }
}
'''
content_make_comp_filename='''
static std::string make_comp_filename(std::string& in_fname) {
  auto dot_pos = in_fname.find_last_of('.');
  if (dot_pos == std::string::npos) {
    return in_fname + "_comp";
  } else {
    return in_fname.substr(0, dot_pos) + "_comp" + in_fname.substr(dot_pos);
  }
}
'''

content_binary_comp='''
static void binary_comp(std::ofstream& fout, std::vector<SegPair>& results) {
  int segN = results.size();
  fout << segN;
  for (auto& result : results) {
    fout << result.first;
    fout << result.second.size();
    for (auto loc : result.second) {
      fout << loc;
    }
  }
  int seq_len = 0;
  std::for_each(bufs_len, bufs_len + threadN,
                [&seq_len](int len) { seq_len += len; });
  fout << seq_len;
  u8 tt;
  int cnt = 0;
  int buf_idx = 0;
  int idx = 0;
  int buf_len = bufs_len[buf_idx];
  char* cursor = bufs[buf_idx];
  while (true) {
    if (cnt == 4) {
      fout << tt;
      cnt = 0;
    }
    if (idx == buf_len) {
      buf_idx++;
      if (buf_idx == threadN) {
        break;
      }
      idx = 0;
      buf_len = bufs_len[buf_idx];
      cursor = bufs[buf_idx];
    }
    tt <<= 2;
    tt += encode_one(*cursor);
    cursor++;
    idx++;
    cnt++;
  }
  fout << tt;
}

'''

content_main1='''
int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: ./comp filename.\\n";
  }
  std::string in_fname = argv[1];
  std::ifstream fin(in_fname);
  if (!fin.is_open()) {
    std::cerr << "[error]: cannot find the file.\\n";
    return 1;
  }
  std::vector<std::thread> threads;
  threadN = -1;
  int tid;
  for (tid = 0; tid < kThreadN && threadN == -1; tid++) {
    threads.push_back(std::thread(thread_routine, tid));
    if (!fin.read(bufs[tid], kMaxBuf)) {
      threadN = tid + 1;
    }
    bufs_len[tid] = fin.gcount();
  }
  fin.close();
  if (threadN == -1) {
    std::cerr << "[error]: sequence two long. please enlarge the buf.\\n";
    return 1;
  }

  MutLock finishN_lck(finishN_mu);
  std::vector<SegPair> results;
  stage = 0;
  for (int i = 0; i < kFindN; i++) {
    finishN = 0;
    MutLock stage_lck(stage_mu);
    stage++;
    stage_cv.notify_all();
    stage_lck.unlock();
    finishN_cv.wait(finishN_lck, []() { return finishN >= threadN; });
    std::cout << "[main]: finish " << i << ".\\n";

    SegMap segmap;
    for (int tid = 0; tid < threadN; tid++) {
      for (auto& result : freqs[tid]) {
        adds_map(segmap, result);
      }
    }
    auto result = std::max_element(segmap.begin(), segmap.end(),
                                   [](const SegPair& a, const SegPair& b) {
                                     return a.second.size() < b.second.size();
                                   });
    freq_result = result->first;
    results.push_back(std::move(*result));
  }
  finishN = 0;
  MutLock stage_lck(stage_mu);
  stage = -1;
  stage_cv.notify_all();
  stage_lck.unlock();
  finishN_cv.wait(finishN_lck, []() { return finishN >= threadN; });

  std::cout << "[main]: finish task!\\n";
  finishN_lck.unlock();
  for (auto& thread : threads) {
    thread.join();
  }
  std::string comp_fname = make_comp_filename(in_fname);
  std::ofstream fout(comp_fname);
  binary_comp(fout, results);
  fout.close();
  return 0;
}
'''
content_include2='''
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
'''
content_Int2String='''
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
	
	str[i] = '\\0';
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
        '''

content_Attach='''
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
'''
content_ReadFile='''
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
            if(str[strlen(str)-1]=='\\n')
                str[strlen(str)-1]='\\0';  //删除\\n
            //存入链表
            int i;
            for(i=0;i<strlen(str);i++)
            {
                if(flag==1)
                {
                    flag=0;
                    prear->num=-1;  //初始化位置
                    Attach('\\n',&prear);
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
'''
content_b_compress='''
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
                case '\\n': i=i;
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
                case '\\n': i=i;
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
'''

content_Findtail='''
void Findtail(Polynomial p,Polynomial *ptail)
{
    int i;
    *ptail=p;
    for(i=0;i<n-1;i++)
    {
        *ptail=(*ptail)->next;
    }
}
'''
content_FindSubstr='''
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
'''
content_DeleteSeed='''
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
'''
content_compare='''
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
'''
content_TentoTwo='''
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
    *(s+k)='\\0';
}
'''
content_CompareStr='''
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
'''
content_find='''
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
        if(pseedl->base=='\\n')
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
    str[i]='\\0';
    free(tmp);

    //格式换位及匹配类型(模糊匹配)
    strcat(str,"\\tApprox\\t");
    
    //存储位置
    char *tmpstr1;
    int digitn;
    tmpstr1=(char *)malloc(sizeof(char)*5);
    Int2String(pstrl->num,tmpstr1,&digitn);  //数字转换为字符串
    
    strcat(str,tmpstr1);  //存储位置
    strcat(str,"\\t");  //格式控制


    //存储长度
    char *tmpstr;
    int digitl;
    tmpstr=(char *)malloc(sizeof(char)*5);
    Int2String(pstrr->num-pstrl->num+1,tmpstr,&digitl);  //数字转换为字符串

    strcat(str,tmpstr);  //存储长度
    strcat(str,"\\t");  //格式控制

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
    str1[i]='\\0';
    free(t);

    Polynomial tt;
    tt=pstrl;
    for(i=0;i<pstrr->num-pstrl->num+1;i++)
    {
        str2[i]=tt->base;
        tt=tt->next;
    }
    str2[i]='\\0';
    free(tt);
    /*比较字符串，存储信息元组*/
    CompareStr(str1,str2,tumple1,tumple2);

    //存储信息元组
    if(tumple1[0]=='e')
    {
        strcat(str,"\\n"); 
    }
    else if(tumple2[0]=='e')
    {
        strcat(str,tumple1); 
        strcat(str,"\\n"); 
    }
    else
    {
        strcat(str,tumple1); 
        strcat(str,"\\t"); 
        strcat(str,tumple2); 
        strcat(str,"\\n"); 
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
'''

content_file_size='''
int file_size(char* filename)
{
    FILE *fp=fopen(filename,"r");
    if(!fp) return -1;
    fseek(fp,0L,SEEK_END);
    int size=ftell(fp);
    fclose(fp);
    
    return size;
}
'''

content_fcompress='''
float fcomprese()
{
    int newnum=0,oldnum=0;  //初始化
    float rate=1.0;  //初始化比例
    newnum=file_size("/home/aubin/bi296/project/dictionary.txt")+file_size("/home/aubin/bi296/project/compress.txt")+file_size("/home/aubin/bi296/project/ids.txt");
    oldnum=file_size("/home/aubin/bi296/project/dna.fasta");
    rate=1.0*newnum/oldnum;
    return rate;
}
'''
content_main2='''
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
    
    
        
        pstr=pseed->next->next->next->next;
        key=FindSubstr(&pseed,&pstr);  //寻找种子及子串
        find(pd,pseed,&pstr);  //延伸并存入离线字典文件
       
        key=FindSubstr(&pseed,&pstr);  //寻找种子及子串
        printf("%d\\n",key);
        find(pd,pseed,&pstr);  //延伸并存入离线字典文件

    
        
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

    printf("%f\\n",frate);
    return 0;
}
'''
content3='''
        #include <stdio>
        int main()
        {
            printf("hello,world");
            return 0;
        }
        '''

intr1='''
        ## SeqCompress

### 一、概述

​        DNA序列压缩算法分为无损和有损两类，可能是垂直、水平或标准模式。SeqCompress算法是一种无损压缩算法，使用算术编码和统计模型来压缩序列，可以应对生物序列的空间复杂性。

### 二、算术编码

​        首先，将输入序列分为两段。一段由FASTA文件标题组成，作为输入序列的标题和标题位置信息。另一段由序列字符组成。

​        由于第二段序列可能包含非ATGC或小写字母，我们需要进一步细分子段，并进行规范化处理。保存该段序列中所有小写字符块的始末位置，并将所有小写字符转换为大写字符。对于非ATGC字符，将这些字符与位置保存在单独的文件中，并从该段序列中删除。

### 三、统计模型

​      （二）中产生的序列作为统计模型的输入。选择长度为l的序列作为输入序列，从中选择长度为n的子段，并保存为单独的文件。该统计模型使用公式（1）计算输入序列中子段s的频率，重复m次，选出m个频率最高的子段。计算压缩率PCR、段压缩率PCR~s~、二进制压缩率PCR~b~ 。如果该段PCR~s~优于PCR~b~，则使用基于段的压缩，否则选择两位二进制压缩，每个子段重复该过程。子段随后从输入序列中删除，剩余的碱基写入两位二进制。过程中生成的所有文件用7z压缩。

公式（1）：
$$
f=\sum_{i=0}^lsi             
$$
公式（2）
$$
PCR=(size of data after compression/size of data before compression)*100
$$
公式（3）
$$
PCR_s=((f*8+n*8)/(f*n*8))*100
$$
公式（4）
$$
PCR_b=((f*n*2)/(f*n*8))*100PCR_b=((f*n*2)/(f*n*8))*100
$$


## 
        '''    
intr2='''
        [streamlit插入网址链接- CSDN搜索](https://so.csdn.net/so/search?q=streamlit插入网址链接&t=&u=)
        该算法使用**统计模型**和**算术编码**。

        **统计模型**找到给定序列中的一个片段的频率，并决定是使用基于片段的压缩还是二进制压缩

        **算术编码**用于编码统计模型中发现的各种片段和剩余核苷酸碱基
        '''  

intr3='''
        
        该算法使用**统计模型**和**算术编码**。

        **统计模型**找到给定序列中的一个片段的频率，并决定是使用基于片段的压缩还是二进制压缩

        **算术编码**用于编码统计模型中发现的各种片段和剩余核苷酸碱基
        '''      


df = pd.DataFrame({
  'alg1': ['SeqCompress','introduction','Pseudocode', 'algorithm'],
  'alg2': ['SeedBaseComp','introduction','Pseudocode', 'algorithm'],
  'alg3': ['Huffman','introduction','Pseudocode', 'algorithm'],
  'result': [' ','result',' ',' '],
  'future': [' ','.',' ',' '],
  'help': [' ','.',' ',' ']
})


line_chart='''
if st.checkbox('line chart'):
    chart_data = pd.DataFrame(
       np.random.randn(20, 3),
       columns=['a', 'b', 'c'])

    st.line_chart(chart_data)
'''

option1 = st.sidebar.selectbox(
    '''
    alg1
    ''',
    df['alg1'].unique()
)

option2 = st.sidebar.selectbox(
    '''
    alg2
    ''',
    df['alg2'].unique()
)

option3 = st.sidebar.selectbox(
    '''
    alg3
    ''',
    df['alg3'].unique()
)

optionres = st.sidebar.selectbox(
    '''
    result
    ''',
    df['result'].unique()
)

option4 = st.sidebar.selectbox(
    '''
    future
    ''',
    df['future'].unique()
)

option5 = st.sidebar.selectbox(
    '''
    help
    ''',
    df['help'].unique()
)

if option1=='SeqCompress' and option2=='SeedBaseComp' and option3=='Huffman' and optionres==' ' and option4==' ' and option5==' ':
    if st.button('Visualization'):
        st.video(video)
    st.write(bk)

if option1=='algorithm':
    if st.checkbox('include'):
        st.code(content_include1,language='c')
    if st.checkbox('encode_one'):
        st.code(content_encode_one,language='c')
    if st.checkbox('encode_forward'):
        st.code(content_encode_forward,language='c')
    if st.checkbox('encode'):
        st.code(content_encode,language='c')
    if st.checkbox('adds_map'):
        st.code(content_adds_map,language='c')
    if st.checkbox('add_map'):
        st.code(content_add_map,language='c')
    if st.checkbox('struct'):
        st.code(content_struct,language='c')
    if st.checkbox('find_most_freq'):
        st.code(content_find_most_freq,language='c')
    if st.checkbox('remove_freq'):
        st.code(content_remove_freq,language='c')
    if st.checkbox('thread_routine'):
        st.code(content_thread_routine,language='c')
    if st.checkbox('make_comp_filename'):
        st.code(content_make_comp_filename,language='c')
    if st.checkbox('binary_comp'):
        st.code(content_binary_comp,language='c')
    if st.checkbox('main'):
        st.code(content_main1,language='c')
if option2=='algorithm':
    if st.checkbox('include'):
        st.code(content_include2,language='c')
    if st.checkbox('Int2String'):
        st.code(content_Int2String,language='c')
    if st.checkbox('Attach'):
        st.code(content_Attach,language='c')
    if st.checkbox('ReadFile'):
        st.code(content_ReadFile,language='c')
    if st.checkbox('b_compress'):
        st.code(content_b_compress,language='c')
    if st.checkbox('add_map'):
        st.code(content_add_map,language='c')
    if st.checkbox('Findtail'):
        st.code(content_Findtail,language='c')
    if st.checkbox('FindSubstr'):
        st.code(content_FindSubstr,language='c')
    if st.checkbox('DeleteSeed'):
        st.code(content_DeleteSeed,language='c')
    if st.checkbox('compare'):
        st.code(content_compare,language='c')
    if st.checkbox('TentoTwo'):
        st.code(content_TentoTwo,language='c')
    if st.checkbox('CompareStr'):
        st.code(content_CompareStr,language='c')
    if st.checkbox('find'):
        st.code(content_find,language='c')
    if st.checkbox('find_size'):
        st.code(content_find_size,language='c')
    if st.checkbox('fcompress'):
        st.code(content_fcompress,language='c')
    if st.checkbox('main'):
        st.code(content_main2,language='c')
if option3=='algorithm':
    if st.checkbox('main'):
        st.code(content1,language='c')
    if st.checkbox('func'):
        st.code(content1,language='c')

if option1=='introduction':
    st.write(intr1)
if option2=='introduction':
    st.write(intr2)
if option3=='introduction':
    st.write(intr3)

if option1=='Pseudocode':
    st.image(pse1,use_column_width=True)
if option2=='Pseudocode':
    st.image(pse2,use_column_width=True)
if option3=='Pseudocode':
    st.image(pse3,use_column_width=True)

if optionres=='result':
    '''
    结果表格
    绘图
    '''

if option4=='.':
    '''
    asd
    '''
    st.write(contentf)

if option5=='.':
    '''
    help
    '''