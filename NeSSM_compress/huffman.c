#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <time.h>

#define NIL_PT		0xffff
#define LEFT_IDX	0
#define RIGHT_IDX   1
#define LEFT_MASK	(1<<LEFT_IDX)
#define RIGHT_MASK	(1<<RIGHT_IDX)
#define SLICE_ARRAY_SIZE(slice_arr)  (slice_arr.tail-slice_arr.head+1)
#define GET_VALUE(ch1,ch2)		 (((UINT16)ch1<<8) + (UINT16)ch2)
#define FREE_CODE_LOW_LIMIT		128

#define CH_TYPE_CTR    0	//control char, the ASCII<32
#define CH_TYPE_SYM	   1	//space, ",","+","(" ":" and so on
#define CH_TYPE_NUM    2	//'0' to '9'
#define CH_TYPE_LET	   3	//letter of the alphabet

#define MAX_LINE 10000

typedef unsigned int  UINT32;
typedef unsigned short UINT16;
typedef unsigned char  BYTE;

typedef struct{
	UINT32 weight;
	char*  huffman_code;
	BYTE   bytes[2];
	UINT16 parent;
	UINT16 l_child;
	UINT16 r_child;
}HM_ENCODE_NODE;

typedef struct{
	UINT32 weight;
	BYTE   bytes[2];
	BYTE   no;
}DOUBLE_BYTES_CODE;

typedef struct{
	DOUBLE_BYTES_CODE *data;
	int len;
}SEARCH_ARRAY;	//double bytes to single byte convert table

typedef struct{
	SEARCH_ARRAY arr;
	int  coding_space;		//2 byte encoding space
	int  encode_mode;		//1:普通方式编码，2：对双字节组合使用单字节编码
	BYTE reserve_char;		//保留字符
}MAP_2BYTES_1BYTE;

typedef struct{
	BYTE flag;
	BYTE left;
	BYTE right;
}HM_DECODE_NODE;

HM_ENCODE_NODE  g_HTree[256*2];			// the global varible, huffman tree used in encoding
char*			g_HM_code_buff=NULL;	// the buffer of "01"string, global varible, be alloc during run time and release when encoding is completed
const BYTE		g_maskArray[]={1,2,4,8,16,32,64,128};
MAP_2BYTES_1BYTE g_transStruct;

int split(char **arr,char *str,char *del)
{
	int i=0;
	char *s=NULL;
	s=strtok(str,del);
	while(s!=NULL)
	{
		*arr++=s;
		s=strtok(NULL,del);
		i++;
	}
	return i;
}

int get_char_type(BYTE ch)
{
	if (ch< ' ')
		return CH_TYPE_CTR;
	else if (ch>='0' && ch<='9')
		return CH_TYPE_NUM;
	else if (ch>='a' && ch<='z'|| ch>='A' && ch<='Z' )
		return CH_TYPE_LET;
	else
		return CH_TYPE_SYM;
}

void init_transStruct()
{
	g_transStruct.arr.data  =NULL;
	g_transStruct.arr.len =0;
	g_transStruct.encode_mode=1;				//1:普通方式编码，2：对双字节组合使用单字节编码
	g_transStruct.reserve_char=0;				//一个特别的字符
}

void free_transStruct()
{
	if ( g_transStruct.arr.data!=NULL)
	{
		free(g_transStruct.arr.data); g_transStruct.arr.data=NULL;
	}

	g_transStruct.arr.len=0;					//the array length
}

/* The compare function be used in qsort*/
int cmp_by_weight(const void *a ,const void *b)
{
	HM_ENCODE_NODE** pp1=(HM_ENCODE_NODE **)a;
	HM_ENCODE_NODE** pp2=(HM_ENCODE_NODE **)b;
	HM_ENCODE_NODE* p1=*pp1;
	HM_ENCODE_NODE* p2=*pp2;

	if ( p1->weight < p2->weight)
		return -1;
	else if ( p1->weight > p2->weight)
		return 1;
	else if ( p1<p2)
		return -1;
	else if (p1>p2)
		return 1;
	else
		return 0;
}

/* The compare function be used in qsort*/
int cmp_by_weight_in_doubleBytesArray(const void *a ,const void *b)
{
	const DOUBLE_BYTES_CODE* p1=a;
	const DOUBLE_BYTES_CODE* p2=b;

	if ( p1->weight < p2->weight)
		return -1;
	else if ( p1->weight > p2->weight)
		return 1;
	else
		return 0;
}

/* The compare function be used in qsort*/
int cmp_by_bytes_in_doubleBytesArray(const void *a ,const void *b)
{
	const DOUBLE_BYTES_CODE* p1=a;
	const DOUBLE_BYTES_CODE* p2=b;
	UINT16 value1=GET_VALUE(p1->bytes[0],p1->bytes[1]);
	UINT16 value2=GET_VALUE(p2->bytes[0],p2->bytes[1]);
	if (value1 < value2)
		return -1;
	else if ( value1 > value2)
		return 1;
	else
		return 0;
}

void dump_Huffman_Tree()
{
	int i;
	printf("no: weight\t ch1 ch2 parent left right huffman_code\n");
	for (i=0;i<512;i++)
	{
		if (g_HTree[i].weight>0)
			printf("%3d:%9d 0x%02x 0x%02x %3d %3d %3d %s\n",i,
				g_HTree[i].weight,g_HTree[i].bytes[0],g_HTree[i].bytes[1],
				g_HTree[i].parent,g_HTree[i].l_child,g_HTree[i].r_child,g_HTree[i].huffman_code);
	}
}


/* print huffman code for every char */
void printf_Huffman_Tree()
{
	int i;
	BYTE ch1,ch2;

	printf("char\tfrequency\thuffman_code\n");
	for (i=0;i<256;i++)
	{
		if (g_HTree[i].weight>0)
		{
			if ( g_HTree[i].bytes[0] == g_transStruct.reserve_char)
			{
				if ( i>= 0x20)
					printf("[\'%c\']",i);
				else
					printf("[%02x]",i);

			}
			else
			{
				ch1=g_HTree[i].bytes[0];
				ch2=g_HTree[i].bytes[1];
				if (ch1 >= ' ' && ch2>= ' ')
				{
					printf("[\'%c%c\']",ch1,ch2);
				}
				else
				{
					printf("[%02x %02x]",ch1,ch2);
				}
			}
			printf("\t%10d\t%s\n",g_HTree[i].weight,g_HTree[i].huffman_code);
		}
	}
}

void search_and_insert( SEARCH_ARRAY *pArray, char ch1,char ch2)
{
	int low,high,mid;
	UINT16  value1,value2,key;
	key=GET_VALUE(ch1,ch2);
	low=0;
	high=pArray->len-1;

	if ( pArray->len==0)
	{
		pArray->data[0].bytes[0]=ch1;
		pArray->data[0].bytes[1]=ch2;
		pArray->data[0].weight=1;
		pArray->len++;
		return;
	}

	value1= GET_VALUE(pArray->data[low].bytes[0], pArray->data[low].bytes[1]);
	value2= GET_VALUE(pArray->data[high].bytes[0],pArray->data[high].bytes[1]);

	if ( key<value1)		//insert the new element before first
	{
		memcpy( pArray->data+1,pArray->data,sizeof(DOUBLE_BYTES_CODE)*(pArray->len));
		pArray->data[low].bytes[0]=ch1;
		pArray->data[low].bytes[1]=ch2;
		pArray->data[low].weight=1;
		pArray->len++;
		return ;
	}
	else if (key>value2)	//insert the new element after last
	{
		pArray->data[high+1].bytes[0]=ch1;
		pArray->data[high+1].bytes[1]=ch2;
		pArray->data[high+1].weight=1;
		pArray->len++;
		return ;
	}

	while (low+1<high)
	{
		mid=(low+high)/2;
		value1= ((UINT16)pArray->data[mid].bytes[0]<<8) + (UINT16)pArray->data[mid].bytes[1];
		if ( key>=value1)
			low=mid;
		else
			high=mid;
	}

	value1= GET_VALUE(pArray->data[low].bytes[0], pArray->data[low].bytes[1]);
	value2= GET_VALUE(pArray->data[high].bytes[0],pArray->data[high].bytes[1]);
	if ( key==value1)
	{
		pArray->data[low].weight++;
		return ;
	}
	else if (key==value2)
	{
		pArray->data[high].weight++;
		return ;
	}
	//insert the new element between low and high
	memcpy( pArray->data+high+1,pArray->data+high,sizeof(DOUBLE_BYTES_CODE)*(pArray->len));
	pArray->data[high].bytes[0]=ch1;
	pArray->data[high].bytes[1]=ch2;
	pArray->data[high].weight=1;
	pArray->len++;
	return ;
}


BYTE getOrder( SEARCH_ARRAY *pArray, char ch1,char ch2,char reserve_char)
{
	int low,high,mid;
	UINT16  value,key;

	low=0;
	high=pArray->len-1;
	key=GET_VALUE(ch1,ch2);
	while (low+1<high)
	{
		mid=(low+high)/2;
		value=GET_VALUE(pArray->data[mid].bytes[0],pArray->data[mid].bytes[1]);
		if ( key>value)
			low=mid+1;
		else if (key<value)
			high=mid-1;
		else
			return pArray->data[mid].no;
	}

	if ( key==GET_VALUE(pArray->data[low].bytes[0], pArray->data[low].bytes[1]))
		return pArray->data[low].no;

	else if (key==GET_VALUE(pArray->data[high].bytes[0], pArray->data[high].bytes[1]))
		return pArray->data[high].no;

	else
		return reserve_char; //means did not find ch1:ch2
}

int packBuffer(char *buff,int len,MAP_2BYTES_1BYTE* pMap)
{
	int i,no,new_len;
	int type1,type2;
	char ch1,ch2;

	for (new_len=0,i=0;i<len-1;)
	{
		ch1=buff[i];
		ch2=buff[i+1];
		type1=get_char_type(ch1);
		type2=get_char_type(ch2);
		no=pMap->reserve_char;
		if (type1==type2)
			no=getOrder( &(pMap->arr),ch1,ch2,pMap->reserve_char);

		if (no !=pMap->reserve_char && type1 == type2)  //replace 2 char with 1 char
		{	buff[new_len++]=no;		i+=2;	}
		else
		{	buff[new_len++]=ch1;	i++;	}
	}
	if (i<len)
		buff[new_len++]= buff[i];
	return new_len;
}

// 返回buff new length
int map_2bytes_to_single_byte(BYTE *buff,int len,MAP_2BYTES_1BYTE *pMap)
{
	int i,j,max_len,new_len;
	int type1,type2;
	void *tmpBuff;

	if ( pMap->arr.data !=NULL)
		free(g_transStruct.arr.data);

	max_len=(256-pMap->coding_space-1)*(256- pMap->coding_space-1);
	pMap->arr.data=(DOUBLE_BYTES_CODE*)malloc(sizeof(DOUBLE_BYTES_CODE)*max_len);

	for (i=0;i<len-1;)
	{
		type1=get_char_type(buff[i]);
		type2=get_char_type(buff[i+1]);
		if (type1==type2)
		{
			search_and_insert( &(pMap->arr),buff[i],buff[i+1]); //make a array with  char pair
			i+=2;
		}
		else
			i++;
	}

	//sort array order by weight
	qsort( pMap->arr.data,pMap->arr.len,sizeof(DOUBLE_BYTES_CODE),cmp_by_weight_in_doubleBytesArray);
	new_len=pMap->coding_space;
	if ( new_len>pMap->arr.len)
		new_len=pMap->arr.len;

	//truncate array g_transStruct.arrWeight
	tmpBuff=malloc(sizeof(DOUBLE_BYTES_CODE)*new_len);
	memcpy(tmpBuff,pMap->arr.data,sizeof(DOUBLE_BYTES_CODE)*new_len);
	free(pMap->arr.data);
	pMap->arr.data=(DOUBLE_BYTES_CODE*)tmpBuff;
	pMap->arr.len=new_len;

	//sort array order by bytes
	qsort( pMap->arr.data,pMap->arr.len,sizeof(DOUBLE_BYTES_CODE),cmp_by_bytes_in_doubleBytesArray);

	for (i=0;i<256;i++)
		g_HTree[i].bytes[0]=pMap->reserve_char; //the all of char is single byte coding

	for (j=0,i=pMap->reserve_char+1;i<256 && j<pMap->arr.len;i++)
	{
		if ( g_HTree[i].weight==0 )
		{
			g_HTree[i].bytes[0]=pMap->arr.data[j].bytes[0];
			g_HTree[i].bytes[1]=pMap->arr.data[j].bytes[1];
			pMap->arr.data[j].no=i;			// i <==> arr.data[j].bytes[0] and bytes[1], single char <==> double char
			j++;
		}
	}

	new_len= packBuffer(buff,len,&g_transStruct);
	for (i=0;i<256;i++)
		g_HTree[i].weight=0;
	for (i=0;i<new_len;i++)
		g_HTree[buff[i]].weight++;
	return new_len;
}

void createHuffmanTree(int mod,char *buff,int len)
{
	//A fix space (256*2 pointer) size slice array
	//The array space is fix, the element count = tail-head+1
	typedef struct
	{
		HM_ENCODE_NODE*  arr[256*2];
		int head;
		int tail;
	}SLICE_ARR;

	SLICE_ARR HT_Arr;
	int new_node_idx,i,j,i1,i2;
	int new_len;	//some of double byte be converted to single char, so the new_len is less than len
	int free_coding_count;	//not used char count,the every free char can denote a double byte combination
	int find_reserve_char;
	char ch;
	BYTE reserve_char;

	memset(g_HTree,0,sizeof(g_HTree));
	for (i=0;i<256;i++)
		g_HTree[i].parent=NIL_PT;

	for (i=0;i<len;i++)
	{
		ch=buff[i];
		g_HTree[ch].weight++;
	}

	find_reserve_char=0;
	free_coding_count=0;
	for (j=0,i=0;i<256;i++)
	{
		if ( g_HTree[i].weight==0 )
		{
			free_coding_count++; //not used char
			if ( !find_reserve_char)
			{
				reserve_char=(BYTE)i;
				find_reserve_char=1;
			}
		}
	}

	for (j=0,i=0;i<256;i++)
	{
		if (g_HTree[i].weight>0)
			HT_Arr.arr[j++]= &(g_HTree[i]); //put the address of element in Huffman Array to array HT_Arr
	}

	HT_Arr.head=0;
	HT_Arr.tail=j-1;
	qsort( HT_Arr.arr,SLICE_ARRAY_SIZE(HT_Arr),sizeof(HM_ENCODE_NODE*),cmp_by_weight);	 //sort array order by weight

	j=257; //the position 256 is reserved for root pointer
	while ( SLICE_ARRAY_SIZE(HT_Arr)>=2)
	{
		i1=HT_Arr.head;
		i2=HT_Arr.head+1;
		if (SLICE_ARRAY_SIZE(HT_Arr)==2)  //the new node is root node
			new_node_idx=256;			//keep the root node is the first place
		else
			new_node_idx=j++;

		g_HTree[new_node_idx].weight=  HT_Arr.arr[i1]->weight +HT_Arr.arr[i2]->weight;
		g_HTree[new_node_idx].l_child= HT_Arr.arr[i1] - g_HTree;
		g_HTree[new_node_idx].r_child= HT_Arr.arr[i2] - g_HTree;
		g_HTree[new_node_idx].parent=NIL_PT;

		HT_Arr.arr[i1]->parent=HT_Arr.arr[i2]->parent=new_node_idx;
		HT_Arr.head+=2;  //remove 2 element from HT_Arr in front

		//insert the new element to proper position
		i=HT_Arr.tail;
		while (i>=HT_Arr.head &&  HT_Arr.arr[i]->weight > g_HTree[new_node_idx].weight)
		{
			HT_Arr.arr[i+1] =HT_Arr.arr[i] ;
			i--;
		}
		HT_Arr.arr[i+1] =  &(g_HTree[new_node_idx]); //Insert this element in postion i+1
		HT_Arr.tail++;
	}
}

void build_huffman_code()
{
	char code[256];
	char *pt;
	int i,j,pos,up;
	int total_char_count=0;
	int total_code_size=0;

	//The first scan
	for (i=0;i<256;i++)
	{
		if (g_HTree[i].weight>0)
		{
			total_char_count++;
			j=i;
			while ( g_HTree[j].parent !=NIL_PT)
			{
				j=g_HTree[j].parent;
				total_code_size++;
			}
		}
	}

	if ( g_HM_code_buff!=NULL)
	{
		free(g_HM_code_buff); g_HM_code_buff=NULL;
	}
	g_HM_code_buff=(char *)malloc(total_char_count + total_code_size); //each of string need a terminal char \0

	//The second scan  encoding
	pt=g_HM_code_buff;
	for (i=0;i<256;i++)
	{
		if (g_HTree[i].weight>0)
		{
			code[sizeof(code)-1]=0;	//the last char be set as 0
			pos=sizeof(code)-2;		//put '0' or '1' to code order by from right to left
			j=i;
			while ( g_HTree[j].parent !=NIL_PT)
			{
				up= g_HTree[j].parent;
				if ( g_HTree[up].l_child==j)
					code[pos--]='0'+LEFT_IDX;
				else
					code[pos--]='0'+RIGHT_IDX;
				j=g_HTree[j].parent;
			}
			strcpy(pt,code+pos+1);
			g_HTree[i].huffman_code=pt;
			pt+= strlen(pt)+1;
		}
	}
}

/*对inFile进行huffman编码，输出到outFile */
void Encoding(int mode,char *inFile,char *outFile)
{
	FILE *fp1=NULL;
	FILE *fp2=NULL;
	UINT32 byteCount;
	char *inBuff=NULL;
	char *outBuff=NULL;
	char head[4];
	if(mode==0)
	{
	    strcpy(head,"IDEN");
	}
	else if(mode==1)
	{
	    strcpy(head,"QVQV");
	}

	fp1=fopen(inFile,"rb");
	if (fp1==NULL)
	{
		printf("can not open file %s\n",inFile);
		goto thisExit;
	}

	fseek(fp1,0,SEEK_END);
	byteCount=ftell(fp1);
	inBuff=(char *)malloc(byteCount);
	fseek(fp1,0,SEEK_SET);
	fread(inBuff,byteCount,1,fp1);
	fclose(fp1); fp1=NULL;

	init_transStruct();
	createHuffmanTree(mode,inBuff,byteCount);
	build_huffman_code();


	{
		UINT32 i,bitIdx;
		UINT32 bitCount=0;
		UINT32 nonLeafNode=0;	//non leaf node count
		UINT32 convert_tab_len=0;	//1 byte code to 2 bytes convert tab length
		char buff[4];

		if (g_transStruct.encode_mode==2)
		{
			i=256-1;
			while ( g_HTree[i].weight==0)
				i--;
			convert_tab_len=i+1;
		}

		for (i=0;i<256*2;i++)
		{
			if ( g_HTree[i].weight>0)
			{
				if (i<256)  //this node is leaf node
					bitCount += (g_HTree[i].weight * strlen(g_HTree[i].huffman_code));
				else        //this node is non-leaf node
					nonLeafNode++;
			}
		}

		//------------------------------------------------------
		fp2=fopen(outFile,"wb");
		if (fp2==NULL)
		{
			printf("can not create file %s\n",outFile);
			goto thisExit;
		}

		//write file header, total 16 bytes
		fwrite(head,4,1,fp2);			//File flag
		fwrite(&nonLeafNode,4,1,fp2);	//Huffman table node count
		fwrite(&byteCount,4,1,fp2);		//The source file byte count
		fwrite(&bitCount,4,1,fp2);		//The bit count of huffman code
		fwrite(&convert_tab_len,4,1,fp2);		//The 1 byte code to 2-bytes convert table length

		memset(buff,0,sizeof(buff));
		buff[0]= g_transStruct.reserve_char;
		fwrite(buff,4,1,fp2);	//The byte count of huffman code

		for (i=0;i<convert_tab_len;i++)
		{
			fwrite(g_HTree[i].bytes,2,1,fp2);	//The byte count of huffman code
		}

		//write huffman table,total nonLeafNode*3 bytes
		for (i=256;i<256*2;i++)
		{
			if ( g_HTree[i].weight>0)
			{
				buff[0]=0;
				if  ( g_HTree[i].l_child <256) //the left node is leaf node
					buff[0] |=LEFT_MASK;

				if  ( g_HTree[i].r_child <256) //the right node is leaf node
					buff[0] |=RIGHT_MASK;

				buff[1]=(BYTE)(g_HTree[i].l_child & 0xff);
				buff[2]=(BYTE)(g_HTree[i].r_child & 0xff);
				fwrite(buff,3,1,fp2);	//The byte count of huffman code
			}
		}

		//write huffman bit data,total (bitCount+7)/8 byte
		outBuff=(char *)malloc((bitCount+7)/8);
		memset(outBuff,0,(bitCount+7)/8);

		i=bitIdx=0;
		while (1)
		{
			BYTE ch;
			BYTE *pHuffmanCode;

			ch=inBuff[i++];
			pHuffmanCode=g_HTree[ch].huffman_code;
			while ( *pHuffmanCode>0)
			{
				if ( *pHuffmanCode++ =='1')
					outBuff[ bitIdx>>3] |= g_maskArray[ bitIdx & 7];
				bitIdx++;
			}
			if (bitIdx>=bitCount)
				break;
		}
		fwrite(outBuff,(bitCount+7)/8,1,fp2);	//The byte count of huffman code
		fclose(fp2); fp2=NULL;					//close file
	}
thisExit:
	free_transStruct();
	if (fp1!=NULL)
	{
		fclose(fp1); fp1=NULL;
	}
	if (inBuff!=NULL)
	{
		free(inBuff); inBuff=NULL;
	}

	if (fp2!=NULL)
	{
		fclose(fp2); fp2=NULL;
	}
	if (outBuff!=NULL)
	{
		free(outBuff); outBuff=NULL;
	}
}

int Write_1ByteOr2Byes( BYTE code, const BYTE *convert_tabl,BYTE reserve_char,char* target)
{
	BYTE ch1,ch2;
	ch1=convert_tabl[code*2];
	ch2=convert_tabl[code*2+1];
	if ( ch1== reserve_char)
	{
		*target=code; return 1;
	}
	else
	{
		target[0]=ch1; target[1]=ch2;	return 2;
	}
}

/*对huffmanFile进行huffman编码，输出到oriFile */
void Decoding(int mode,char *huffmanFile,char *oriFile)
{
	FILE *fp1=NULL;
	FILE *fp2=NULL;
	UINT32 *pDWORD;
	UINT32 i=0,nonLeafNode;	//nonLeafNode, the non leaf node in huffman tree
	UINT32 bitCount;		//the bit count of huffman code
	UINT32 bitIdx;
	UINT32 oriByteCount;	//original file byte count
	UINT32 convert_tab_len;
	BYTE currNode;			//the index of next node or char ASCII
	BYTE reserve_byte;
	BYTE convert_tab[256*2];
	char head[4];

	HM_DECODE_NODE huffmanTree[256];  //huffmanTree, only the non-leaf node be store in this array
	char buff[32];  //changed

	char *huffmanData=NULL;
	char *outBuff=NULL;

    if(mode==0)
    {
        strcpy(head,"IDEN");
    }
    else if(mode==1)
    {
        strcpy(head,"QVQV");
    }
    else printf("error!\n");

	fp1=fopen(huffmanFile,"rb");
	if (fp1==NULL)
	{
		printf("can not open file %s\n",huffmanFile);
		goto thisExit;
	}

	//read file header, total 16 bytes
	fread(buff,24,1,fp1);	//File flag
	if (strncmp(buff,"IDEN",4)==0)
	{
		printf("Decoding the identifier\n");
	}
	else if(strncmp(buff,"QVQV",4)==0)
	{
	    printf("Decoding the qulity value\n");
	}
	else
	{
	    printf("The file %s is not huffman zip file\n",huffmanFile);
		goto thisExit;
	}
	pDWORD= (UINT32*)(buff+4);	nonLeafNode=*pDWORD;
	pDWORD= (UINT32*)(buff+8);	oriByteCount=*pDWORD;
	pDWORD= (UINT32*)(buff+12); bitCount=*pDWORD;
	pDWORD= (UINT32*)(buff+16); convert_tab_len=*pDWORD;
	printf("conver_tab_len:%d\n",convert_tab_len);
	reserve_byte=(BYTE)(buff[20]);

	memset(convert_tab,reserve_byte,sizeof(convert_tab));
	if (convert_tab_len>0)
		fread((void *)convert_tab,convert_tab_len*2,1,fp1);

	huffmanData=(char *)malloc((bitCount+7)/8);

	outBuff=(char *)malloc(oriByteCount);

	for (i=0;i<nonLeafNode;i++)  //changed
	{
		int t=fread(buff,3,1,fp1);	//File flag
		huffmanTree[i].flag =buff[0];
		huffmanTree[i].left =buff[1];
		huffmanTree[i].right =buff[2];
	}
	int t =fread(huffmanData,(bitCount+7)/8,1,fp1);
	fclose(fp1); fp1=NULL;

	for (i=0,currNode=0,bitIdx=0;bitIdx<bitCount;bitIdx++)
	{
		BYTE flag;
		BYTE currBit;

		currBit=huffmanData[ bitIdx>>3] & g_maskArray[ bitIdx & 7];
		flag=huffmanTree[currNode].flag;
		if ( currBit )  //the current bit is 1, goto right branch
		{
			if ( !(flag & RIGHT_MASK))
				currNode=huffmanTree[currNode].right; //
			else
			{
				//In this case, the right field is ASCII of huffman code
				int bytes=Write_1ByteOr2Byes( huffmanTree[currNode].right, convert_tab,reserve_byte,outBuff+i);
				i+=bytes;
				currNode=0;
			}
		}
		else  //the current bit is 0, goto left branch
		{
			if (!(flag & LEFT_MASK))
				currNode=huffmanTree[currNode].left;
			else
			{
				//In this case, the left field is ASCII of huffman code
				int bytes=Write_1ByteOr2Byes( huffmanTree[currNode].left, convert_tab,reserve_byte,outBuff+i);
				i+=bytes;
				currNode=0;
			}
		}
	}

	if (i != oriByteCount)
	{
		printf("Error data\n"); goto thisExit;
	}
	fp2=fopen(oriFile,"wb");//open output file
	if (fp2==NULL)
	{
		printf("can not create file %s\n",oriFile);
		goto thisExit;
	}
	fwrite(outBuff,oriByteCount,1,fp2);
	fclose(fp2); fp2=NULL; //close output file

thisExit:
	if (fp1!=NULL)
	{	fclose(fp1); fp1=NULL; }
	if (huffmanData!=NULL)
	{		free(huffmanData); huffmanData=NULL;}
	if (fp2!=NULL)
	{		fclose(fp2); fp2=NULL; }
	if (outBuff!=NULL)
	{	free(outBuff); outBuff=NULL; }
	if ( g_HM_code_buff!=NULL)
	{	free(outBuff); g_HM_code_buff=NULL; }
}

void minus_strand(const char *buffle,char *minus_buffle,int len)
{
	int circle_f;
	for(circle_f=0;circle_f<len;circle_f++)
	{
		if (buffle[circle_f]=='A')
			minus_buffle[len-1-circle_f]='T';
		else if (buffle[circle_f]=='T')
			minus_buffle[len-1-circle_f]='A';
		else if (buffle[circle_f]=='C')
			minus_buffle[len-1-circle_f]='G';
		else
			minus_buffle[len-1-circle_f]='C';
	}
}

void zipfile(char *infile,char *outfile)
{
    FILE *fin = NULL;
    FILE *fout = NULL;
    int i=0;

    char info[MAX_LINE];
    char sequence[MAX_LINE];
    char plus[4];
    char QValue[MAX_LINE];

    UINT32 record_number=0;
    char *QV_buffer=NULL;
    char *INFO_buffer[10000];
    char *name=NULL;
    char *last_name=NULL;
    int flag=0; //first in a continues
    int count=0;
    char tmp[10000];
    char other[10000];
    char error[10000];

    fin=fopen(infile,"rb");
	if (fin==NULL)
	{
		printf("can not open file %s\n",infile);
		//goto thisExit;
	}

	fout=fopen(outfile,"wb");
    if (fout==NULL)
    {
        printf("can not create file %s\n",outfile);
        //goto thisExit;
    }
    QV_buffer=(char *)malloc(MAX_LINE*100000);
    last_name=(char *)malloc(MAX_LINE);
    memset(last_name,0,MAX_LINE);
    name=(char *)malloc(MAX_LINE);

    record_number=1;

    while(fgets(info,MAX_LINE,fin)!=NULL)
    {

        printf("records: %d\n",record_number);
        strcpy(tmp,info);
		i=split(INFO_buffer,tmp,"|");
        memset(name,0,MAX_LINE);
        memset(other,0,MAX_LINE);
        strcpy(name,INFO_buffer[0]);

		for(count=1;count<5;count++)
		{
		    strcat(name,"|");
		    strcat(name,INFO_buffer[count]);
		}
		for(count=5;count<i;count++)
		{
		    strcat(other,"|");
		    strcat(other,INFO_buffer[count]);
		}
		if(strcmp(last_name,name)==0) // not first record
		{
		    fprintf(fout,"%s",other); //write part one
		}
		else  //not first record
		{
		    fprintf(fout,"%s",info);  //write complete one
		}

		memset(last_name,0,MAX_LINE);
		strcpy(last_name,name);
		memset(info,0,MAX_LINE);
        fgets(sequence,MAX_LINE,fin);
		fgets(plus,4,fin);  //+
		fgets(QValue,MAX_LINE,fin);  //qulity value
        strcat(QV_buffer,QValue);
		record_number++;
    }

    char *outFile2="QValue";
    FILE *fout2;
    fout2=fopen(outFile2,"wb");
    if (fout2==NULL)
    {
        printf("can not create file %s\n",outFile2);
        //goto thisExit;
    }
    fprintf(fout2,"%s",QV_buffer);
    fclose(fout2);
    fclose(fout);
    //char *outFile4="QValue";
    char *outFile3="QV.hm";
    char *outFile4="INFO.hm";
    Encoding(0,outfile,outFile4);
    Encoding(1,outFile2,outFile3);
}

void unzipfile(char *infile,char *outfile,char *index_file)
{
    FILE *fin=NULL,*fout=NULL,*findex=NULL;
    FILE *finfo=NULL,*fqv=NULL;
    //for test
    char *file1="INFO.hm";
    char *file2="QV.hm";
    char *file3="INFO.out";
    char *file4="QV.out";
    //variable for records
    char *INFO_buffer[MAX_LINE];
    char *INFO_buffer2[MAX_LINE];
    char *errors[3];
    char base1,base2;
    char *name,*last_name,*identifier;
    char info[MAX_LINE];
    char tmp[MAX_LINE];
    char qvalue[MAX_LINE];
    int gi_record;
    int position=0;
    char specie[200];
    char line2[MAX_LINE];
    int Flength,Tlength;

    //variable for index_file
    FILE *fpath=NULL;
    char index_line[MAX_LINE];
    char *index_name[MAX_LINE];
    char *path=NULL;
    char *buffle=NULL,*minus_buffle=NULL;
    char *tmp_buffle=NULL;
    char *ref_name=NULL,*sequence_name;

    int len,gi_index,flag,len_buf,point;
    int len_sequence_name;

    int i=0,count=0,a=0,j;
    //for test
    printf("-----------------------------------------------get info----------------------------------------\n");
    Decoding(0,file1,file3);
    Decoding(1,file2,file4);
    printf("--------------------------------------------recover the records----------------------------------------\n");
    //open files that needed
    fout=fopen(outfile,"w");
    if(fout==NULL){printf("Warning:can not open %s file in this directory.\n",outfile);exit(1);}
    path=(char *)malloc(2000*sizeof(char));
    finfo=fopen(file3,"r");
    if(finfo==NULL){printf("Warning:can not find the info file.\n");exit(1);}
    fqv=fopen(file4,"r");
    if(fqv==NULL){printf("Warning:can not find the qulity_value file.\n");exit(1);}
    findex=fopen(index_file,"r");
    if(findex==NULL){printf("Warning:can not find the index file.\n");exit(1);}

    //apply mem for recover
    int records=0;
    name=(char *)malloc(MAX_LINE);
    last_name=(char *)malloc(MAX_LINE);
    identifier=(char *)malloc(MAX_LINE);
    //printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ready!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    while(fgets(info,MAX_LINE,finfo)!=NULL)
    {
        records++;
        printf("records %d\n",records);
        //printf("TEST 1\n");
        memset(qvalue,'\0',MAX_LINE);
        memset(line2,'\0',MAX_LINE);
        fgets(qvalue,MAX_LINE,fqv);
        Tlength=strlen(qvalue)-1;
        //printf("-----------------------------------------------QV length is :%d--------------------------\n",Tlength);
        //printf("TEST 2\n");
        if(strncmp(info,"@",1)==0)  //new record series & complete identifier
        {
            //printf("TEST 3\n");
            strcpy(identifier,info);
            i=split(INFO_buffer,info,"|");
           // printf("TEST \n");
            strcpy(name,INFO_buffer[0]);
            gi_record=atoi(INFO_buffer[1]);
           // printf("TEST 4\n");
            memset(specie,'\0',200);
            strcpy(specie,INFO_buffer[4]);
            for(count=1;count<5;count++)  //name
            {
                strcat(name,"|");
                strcat(name,INFO_buffer[count]);
            }
            printf("%s\n",name);
            flag=0;
            memset(index_line,'\0',MAX_LINE);
            while(fgets(index_line,MAX_LINE,findex)!=NULL)					//get the specie's information in the index
            {
                //printf("TEST 5\n");
                index_line[strlen(index_line)-1]='\0';
                split(index_name,index_line,"\t");
                if(memcmp(specie,index_name[3],strlen(specie))==0||(atoi(index_name[0])==gi_record))  //add one condition
                {
                    flag=1;
                    len=atoi(index_name[1]);
                    memset(path,'\0',2000*sizeof(char));
                    strcat(path,index_name[4]);
                    gi_index=atoi(index_name[0]);
                    break;
                }
            }
            rewind(findex);
            if(flag==0)
            {
                printf("Warning:%s can't find in the index!!!",specie);
                continue;
            }

            //recover the record;
            //printf("%s\n",path);
            fpath=fopen(path,"r");
            if(fpath==NULL)
            {
                printf("Warning:can not find the index file pathway in this directory.\n");
                exit(1);
            }
            buffle=(char *)malloc(len*sizeof(char));
            memset(buffle,'\0',len*sizeof(char));

            tmp_buffle=(char *)malloc(len*2*sizeof(char));
            ref_name=(char *)malloc(400*sizeof(char));
            memset(tmp_buffle,'\0',len*2*sizeof(char));
            memset(ref_name,'\0',400*sizeof(char));
            fread(tmp_buffle,len*2,1,fpath);
            fclose(fpath);
            flag=0;
            char *p1,*p2,*p3;
            p1=tmp_buffle;p2=ref_name;p3=buffle;

            while(*p1 != '\0')			//get the specie's name and its sequence
            {
                if(*p1 == '\n') {p1++;flag+=1;continue;}
                else
                {
                    if(flag==0){*p2++=*p1++;}
                    else{*p3++=*p1++;}
                }
            }
            free(tmp_buffle);
            len_buf=strlen(ref_name);
            sequence_name=(char *)malloc((len_buf+2)*sizeof(char));
            memset(sequence_name,'\0',(len_buf+2)*sizeof(char));
            sequence_name[0]='@';
            strcat(sequence_name,ref_name);
            len_sequence_name=strlen(sequence_name);
            minus_buffle=(char *)malloc(len*sizeof(char));
            position=atoi(INFO_buffer[5]);
            minus_strand(buffle,minus_buffle,len);


            //do the length
            Flength=Tlength;

            for(count=7;count<i;count++)
            {
                //strcpy(INFO_buffer2[count],INFO_buffer[count]);
                INFO_buffer2[count]=INFO_buffer[count];
                //printf("info buffer2%s\n",INFO_buffer2[count]);
            //printf("1111 info buffer  %s\n",INFO_buffer[count]);
                split(errors,INFO_buffer2[count],":");
                if(strncmp(errors[1],"-",1)==0) //insertion
                {
                    Flength--;
                    //printf("insertion\t");
                }
                else if(strncmp(errors[2],"-",1)==0) //deletion
                {
                    Flength++;
                    //printf("deletion\t");
                }
            }
            //printf("\n");

            //printf("2222 info buffer  %s\n",INFO_buffer[7]);

            //printf("-------------------False length is :%d-----------------------\n",Flength);
            if(strncmp(INFO_buffer[6],"+",sizeof(char))==0)
            {
                for(a=0;a<Flength;a++)
                {
                    line2[a]=buffle[a+position];
                }
            }
            else if(strncmp(INFO_buffer[6],"-",sizeof(char))==0)
            {
                for(a=0;a<Flength;a++)
                {
                    line2[a]=minus_buffle[a+(len-position-1)];
                }
            }

            //add errors
            for(count=7;count<i;count++)
            {
                //printf("info %s\n",INFO_buffer[count]);
                split(errors,INFO_buffer[count],":");
                //printf("info after num%s\n",INFO_buffer[count]+strlen(INFO_buffer[count])+1);
                point=atoi(errors[0]);
                errors[1] = INFO_buffer[count]+strlen(INFO_buffer[count])+1;
                base1=errors[1][0];
                errors[2] = errors[1]+2;
                base2=errors[2][0];
                //printf("%s\t%s\t%s\n",errors[0],errors[1],errors[2]);
                //printf("%s\t%s\t%s\n",errors[0],errors[1],errors[2]);
                //printf("%s\t%s\t%s\n",errors[0],base1,base2);
                if(strncmp(errors[2],"-",1)==0) //deletion
                {
                    for(j=point;j<strlen(line2);j++)
                    {
                        line2[j]=line2[j+1];
                    }
                    line2[strlen(line2)]='\0';
                    //printf("deletion %d\n",strlen(line2));
                }
                else if(strncmp(errors[1],"-",1)==0) //insertion
                {
                    for(j=strlen(line2)+1;j>point;j--)
                    {
                        line2[j]=line2[j-1];
                    }
                    line2[point]=base2;
                    //printf("insertion %d\n",strlen(line2));
                }
                else //substitution
                {
                    line2[point]=base2;
                    //printf("substitution %d\n",strlen(line2));
                }
            }
            //printf("--------------------------------------------Ture length is :%d----------------------------\n",strlen(line2));

            fprintf(fout,"%s",identifier);
            fprintf(fout,"%s\n",line2);
            fprintf(fout,"+\n");
            fprintf(fout,"%s",qvalue);
        }
        else                          //part identifier
        {
            memset(identifier,'\0',MAX_LINE);
            strcpy(identifier,name);
            //strcat(identifier,"|");
            strcat(identifier,info);
            i=split(INFO_buffer,info,"|");
            position=atoi(INFO_buffer[0]);
            //strcpy(strand,INFO_buffer[1]);

            //do the length
            Flength=Tlength;
            for(count=2;count<i;count++)
            {
                //strcpy(INFO_buffer2[count],INFO_buffer[count]);
                INFO_buffer2[count]=INFO_buffer[count];
                split(errors,INFO_buffer2[count],":");
                if(strncmp(errors[1],"-",1)==0) //insertion
                {
                    Flength--;
                    //printf("insertion\t");
                }
                else if(strncmp(errors[2],"-",1)==0) //deletion
                {
                    Flength++;
                    //printf("deletion\t");
                }
            }
            //printf("-------------------False length is :%d-------------------\n",Flength);

            if(strncmp(INFO_buffer[1],"+",sizeof(char))==0)
            {
                for(a=0;a<Flength;a++)
                {
                    line2[a]=buffle[a+position];
                }
            }
            else if(strncmp(INFO_buffer[1],"-",sizeof(char))==0)
            {
                for(a=0;a<Flength;a++)
                {
                    line2[a]=minus_buffle[a+(len-position-1)];
                }
            }

            //add error

            for(count=2;count<i;count++)
            {
                split(errors,INFO_buffer[count],":");
                point=atoi(errors[0]);
                errors[1] = INFO_buffer[count]+strlen(INFO_buffer[count])+1;
                base1=errors[1][0];
                errors[2] = errors[1]+2;
                base2=errors[2][0];
                //printf("%s\t%s\t%s\n",errors[0],errors[1],errors[2]);
                if(strncmp(errors[2],"-",1)==0) //deletion
                {
                    printf("%s:%s:%s\n",errors[0],errors[1],errors[2]);
                    for(j=point;j<strlen(line2);j++)
                    {
                        line2[j]=line2[j+1];
                    }
                    line2[strlen(line2)]='\0';
                    //printf("deletion %d\n",strlen(line2));
                }
                else if(strncmp(errors[1],"-",1)==0) //insertion
                {
                    for(j=strlen(line2)+1;j>point;j--)
                    {
                        line2[j]=line2[j-1];
                    }
                    line2[point]=base2;
                    //printf("insertion %d\n",strlen(line2));
                }
                else //substitution
                {
                    line2[point]=base2;
                    //printf("substitution %d\n",strlen(line2));
                }
            }
            //printf("----------------------------------------------Ture length is :%d------------\n",strlen(line2));
            fprintf(fout,"%s",identifier);
            fprintf(fout,"%s\n",line2);
            fprintf(fout,"+\n");
            fprintf(fout,"%s",qvalue);
        }
        memset(info,'\0',MAX_LINE);
        memset(qvalue,'\0',MAX_LINE);
    }
    /*
    remove("QV.hm");
    remove("QV.out");
    remove("QValue");
    remove("INFO.hm");
    remove("INFO.out");
    */
}


int main(int argc, char* argv[])
{
	int mode=1;
	int bEncoding=1;
	int i;
	char zipFile[260];
	char srcFile[260];
	char refFile[260];

	if ( argc<4)
	{
		printf("\t%s a zipfile sourcefile\nor\n\t%s e zipfile sourcefile indexfile\n",argv[0],argv[0]);
		return 1;
	}

	if  ( strcmp(argv[1],"a")==0)  //encoding
		bEncoding=1;
	else if  (strcmp(argv[1],"e")==0)  //decoding
		bEncoding=0;
	else
	{
		printf("The first command must be a or e\n");
		return 1;
	}

	zipFile[0]=0;
	srcFile[0]=0;
	refFile[0]=0;

	for (i=2;i<argc;i++)
	{
		if(bEncoding)
		{
		    if (strlen(zipFile)==0)
				strcpy(zipFile,argv[i]);

			else if (strlen(srcFile)==0)
				strcpy(srcFile,argv[i]);
		}
		else
		{
		    if (strlen(zipFile)==0)
				strcpy(zipFile,argv[i]);

			else if (strlen(srcFile)==0)
				strcpy(srcFile,argv[i]);
            else if (strlen(refFile)==0)
                strcpy(refFile,argv[i]);
		}

	}

	if ( bEncoding==1 && (strlen(zipFile)==0 || strlen(srcFile)==0))
	{
		printf("Must enter The zipfile name and source file\n");
		return 1;
	}
	else if(bEncoding==0 && (strlen(zipFile)==0 || strlen(srcFile)==0 || strlen(refFile)==0))
	{
	    printf("Must enter The zipfile name, source file and index file\n");
		return 1;
	}

    time_t start,end;
	start=time(NULL);
    printf("Start...\n");

	if ( bEncoding)
		zipfile(srcFile,zipFile);
	else
		unzipfile(zipFile,srcFile,refFile);

    printf("The end!\n");
    end=time(NULL);
    printf("print %f s\n",difftime(end,start));
    return 0;
}
