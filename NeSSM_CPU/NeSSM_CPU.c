#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define MAX_LINE 10000
#define NN 200000
#define PI 3.1415926
//======================================================================================
//function:usage
//======================================================================================
void usage()
{
	printf("This is a customizable metagenome simulation system: NeSSM (Next-generation Sequencing Simulator for Metagenomics).\n");
	printf("By combining completed genomes currently available and some initial metagenome sequencing information, it can create a simulated metagenome closer-than-ever to the reality.\n");
	printf("With this system, an optimized metagenomics analysis strategy for a specific project can be created according to the sequencing platforms available to the user, the complexity of tared metagenome, and the goal of metagenome sequencing. Currently, NeSSM supports 454,sange and Illumina sequencing platforms.\n");
	printf("NeSSM can be helpful to evaluate or develop tools for some metagenomics applications.\n\n");
	printf("Usage: ./NeSSM_CPU  [options]\n");
	printf("[options]:\n");
	printf("-t <simulation type>            'genome' or '16s'. simulate genome sequence or 16s rRNA sequence, default genome sequence.\n");
	printf("-r <reads_number>               number of reads to generate, default 1000.\n");
	printf("-o <output_file>                output file.\n");
	printf("-m <sequencing_method>          'sanger' or '454' or 'illumina' or 'pacbio', default 'illumina'.\n");
	printf("-e <single or paired-end>       'single'(single-end) or 'pair'(paired-end), for both 454 and illumina sequencing,\n");
	printf("                                0 means 'single' and 1 means 'pair'. default 0.\n");
	printf("-list <input_file>              16s sequence file or a metagenome composition table.\n");
	printf("-index <index_file>             generate a set of NGS reads for <index_file>. Pay attantion to the genome index or the 16S index.\n" );
	printf("-l <sequencing length>          length of the reads, default for illumina sequencing 36(36bp) for pacbio is 3000bp.\n");
	printf("-w <width of pairs>             gap length of pairs, for paired-end only. default 200.\n");
	printf("-c <error_model_file>           error model file.default 'simulation2.config'.\n");
	printf("-exact <exact simulation>       0 means the length of read is decided by the parameter -l, 1 means the length of read is decided by the distribution of length according to a real data, default is 0.\n");
	printf("-buff <number_buff>             the buff number in memory.default 1.\n" );
	printf("-b <coverage_bias_file>         the coverage bias file\n");
}
void runNeSSM(int argc,char **argv);   //declare a function:runNeSSM
//======================================================================================
//function:rand_length--generate a read's length according to average and standard deviation,
//		   the result obeys normal distribution. the length doesn't exceed the max.
//=====================================================================================
int rand_length(float average,float sd,int max,int min)
{
	float out,random,test;
	int flag=0;

	if(sd==0)
	{
		return (int)average;
	}

	while(flag==0)
	{
		random=(float)((float)rand()/RAND_MAX-0.5)*50*sd+average;
		test=(float)rand()/RAND_MAX*2.0/sqrt(2.0*3.1415926)/sd;
		out=random;
		random=1.0/sqrt(2.0*3.1415926)/sd*(float)exp((-1)*(random-average)*(random-average)/2.0/sd/sd);

		if((test<random)&&(out>min)&&(out<=max))    //round-off
		{
			flag=1;
		}
	}
	return (int)(out);
}

//======================================================================================
//function:rand_length2--generate a read's length according to average and standard deviation,
//		   the result obeys normal distribution. For pacbio.
//=====================================================================================
double AverageRandom(double min,double max)//产生(min,max)之间均匀分布的随机数
{
        int MINnteger = (int)(min*10000);
        int MAXnteger = (int)(max*10000);
        int randInteger = rand()*rand();
        int diffInteger = MAXnteger - MINnteger;
        int resultInteger = randInteger % diffInteger + MINnteger;
        return resultInteger/10000.0;
}
double Normal(double U1,double U2)
{
        double Z = sqrt(-2*log(U1)) * cos(2*PI*U2);
        return Z;
}
int rand_length2(int length,char use,int max_length) //normal distribution
{
	    double X;
        double U1,U2;
        double Z;
        double SIGMA=log(length);
        double MIU;
        if(use == '1'){MIU=0.22;} //for pair-end gap length
        else if (use == '0')    //for pacbio reads length
        {
            MIU=0.8;
        }
        do{
        U1 = AverageRandom(0,1);
        U2 = AverageRandom(0,1);
        Z = Normal(U1,U2);
        X = SIGMA + (Z * MIU);
        }while((int)(exp(X)) <= 0 || (int)(exp(X)) >= max_length);
        return (int)(exp(X));
}

//======================================================================================
//function:split--split a string according with a small string and return how many segments
//          after spliting.
//======================================================================================
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
//======================================================================================
//function:sub_single--substitute one base
//======================================================================================
char sub_single(char b1,char b2,char b3,float p1,float p2,float p3)
{
	float per1,per2,per;
	per=(float)(rand())/(float)(RAND_MAX);
	per1=p1/(p1+p2+p3);
	per2=p2/(p1+p2+p3);
    if(per<=per1)
		return b1;
	else if(per>per1 && per<=(per1+per2))
		return b2;
	else
		return b3;
}
//======================================================================================
//function:sub_all--substitute the inputed base for DNA
//======================================================================================
char sub_all(char seq,float ac,float ag,float at,float cg,float ct,float ca,float gt,float ga,float gc,float ta,float tc,float tg)
{
    if(seq=='A')
		return sub_single('C','G','T',ac,ag,at);
	else if(seq=='C')
		return sub_single('G','T','A',cg,ct,ca);
	else if(seq=='G')
		return sub_single('T','A','C',gt,ga,gc);
	else if(seq=='T')
		return sub_single('A','C','G',ta,tc,tg);
}
//======================================================================================
//function:sub_all_r--substitute the inputed base for RNA
//======================================================================================
char sub_all_r(char seq,float ac,float ag,float at,float cg,float ct,float ca,float gt,float ga,float gc,float ta,float tc,float tg)
{
    if(seq=='A')
		return sub_single('C','G','T',ac,ag,at);
	else if(seq=='C')
		return sub_single('G','T','A',cg,ct,ca);
	else if(seq=='G')
		return sub_single('T','A','C',gt,ga,gc);
	else if(seq=='T')
		return sub_single('A','C','G',ta,tc,tg);
}
//======================================================================================
//function:ins_all--generate a base for insertion for DNA
//======================================================================================
char ins_all()
{
	float random=(float)(rand())/(float)(RAND_MAX);
    if(random<=0.25)
		return 'A';
	else if(random<=0.5)
		return 'C';
	else if(random<=0.75)
		return 'G';
	else
		return 'T';
}
//======================================================================================
//function:ins_all_r--generate a base for insertion for RNA
//======================================================================================
char ins_all_r()
{
	float random=(float)(rand())/(float)(RAND_MAX);
    if(random<=0.25)
		return 'A';
	else if(random<=0.5)
		return 'C';
	else if(random<=0.75)
		return 'G';
	else
		return 'T';
}
//======================================================================================
//main function
//======================================================================================
main(int argc,char **argv)
{
	if(argc>=6)
	{
		int t,lable;
		lable=0;
		for(t=1;t<=argc-1;t+=2)
		{
			if(strcmp(argv[t],"-index")==0){lable+=1;}
			else if(strcmp(argv[t],"-list")==0){lable+=1;}
			else if(strcmp(argv[t],"-o")==0){lable+=1;}
		}
		if(lable==3)
		{
			time_t start,end;
			start=time(NULL);
			printf("This is a CPU program!\n");
			printf("Start...\n");
			runNeSSM(argc,argv);
			printf("The end!\n");
			end=time(NULL);
			printf("print %f s\n",difftime(end,start));
			return 0;
		}
		else
		{
			printf("You must set the parameter '-o','list','index'.\n\n");
			usage();
			return -1;
		}
	}
	else
	{
		usage();
		return -1;
	}
 }
//======================================================================================
//function:runNeSSM--the main function to simulate
//======================================================================================
void runNeSSM(int argc,char **argv)
{
	//sub-function
	int bp();
	void re_everyread();
	void read_print();
	void input_right();
	void minus_strand();
	void minus_strand_r();
	int species_number();
	int get_qv_type();
	int first_line();
	int remain_lines();
	int quality();
	int length_distribution();
	int start_position();
	int count_total_copy();
	void checkposition();

	//some variants
	int reads,gap,length,circle,lable,local,flag_buff,buff,gap_lable,lable_method,len,len_sequence_name,everynumber[2000],copy_number[2000],min,exact,flag_bias,copy,totalcopy;
	int species_reads,len_buf,i,x,y,flag,total,point1,point2,point,flag_pn,buff_point,buff_point1,buff_point2,length_max,every_length;
	buff=1;
	gap_lable=0;
	total=0;
	reads=1000;
	gap=200;
	length=50;
	lable=0;
	exact=0;
	char index[MAX_LINE],list[MAX_LINE],fastq[MAX_LINE],name[MAX_LINE],*index_name[MAX_LINE],*path,*ref_name,*number[2],*tmp_buffle;
	char *buffle,*sequence_name,name1[2000][200],*line2_1,*line2_2,*line4_1,*line4_2,*line2,*line4,*minus_buffle,bias[MAX_LINE];
	char config[MAX_LINE]="simulation.config";
	char method[20]="illumina";
	char simu_type[20]="gemone";
	FILE *fp,*findex,*fpath,*fs1,*fs2,*fs,*fc,*fbias;
	float random_point,reads_number,every_percent,sd_ratio[1],type[3],sd,subratio[12],*err_rand,length_rand[4000];
	int gi,max_t,get_gi,flag_gi;
	char *all,*part_t[2],*value_t[NN];
	float value_f[NN],qual[50];

	flag_bias=0;   //don't use the coverage bias file

	for(i=0;i<50;i++)
		qual[i]=pow(10,(i/(-10.0)));

	for(circle=1;circle<=argc-1;circle+=2)					//read the input commend
	{
		if(strcmp(argv[circle],"-r")==0){reads=atoi(argv[circle+1]);}
		else if(strcmp(argv[circle],"-index")==0){strcpy(index,argv[circle+1]);}
		else if(strcmp(argv[circle],"-e")==0){lable=atoi(argv[circle+1]);}
		else if(strcmp(argv[circle],"-list")==0){strcpy(list,argv[circle+1]);}
		else if(strcmp(argv[circle],"-m")==0){strcpy(method,argv[circle+1]);if(strstr(method,"pacbio")!=NULL){length=9000;}}
		else if(strcmp(argv[circle],"-l")==0){length=atoi(argv[circle+1]);}
		else if(strcmp(argv[circle],"-o")==0){local=circle;}
		else if(strcmp(argv[circle],"-w")==0){gap=atoi(argv[circle+1]);gap_lable=1;}
		else if(strcmp(argv[circle],"-c")==0){strcpy(config,argv[circle+1]);}
		else if(strcmp(argv[circle],"-buff")==0){buff=atoi(argv[circle+1]);}
		else if(strcmp(argv[circle],"-exact")==0){exact=atoi(argv[circle+1]);}
		else if(strcmp(argv[circle],"-b")==0){strcpy(bias,argv[circle+1]);flag_bias=1;}
		else if(strcmp(argv[circle],"-t")==0){strcpy(simu_type,argv[circle+1]);}
		else {printf("Waring:there is no this parameter '%s',please see the Usage again!\n",argv[circle]);exit(1);}
	}

	if(lable==0)	//single							//open the file to reserve the simulation
	{
		memset(fastq,'\0',MAX_LINE*sizeof(char));
        strcpy(fastq,argv[local+1]);
        strcat(fastq,".fq");
		fs=fopen(fastq,"w");
		if(fs==NULL){printf("Warning:can not make the fastq file in this directory.\n");exit(1);}
		printf("the fastq file:%s.\n",fastq);
	}
	if(lable==1)      //pair-end
	{
		memset(fastq,'\0',MAX_LINE*sizeof(char));
        strcpy(fastq,argv[local+1]);
		strcat(fastq,"_1.fq");
		fs1=fopen(fastq,"w");
		if(fs1==NULL){printf("Warning:can not make the fastq_1 file in this directory.\n");exit(1);}
		printf("the fastq_1 file:%s.\n",fastq);

		memset(fastq,'\0',MAX_LINE*sizeof(char));
		strcpy(fastq,argv[local+1]);
		strcat(fastq,"_2.fq");
        fs2=fopen(fastq,"w");
		if(fs2==NULL){printf("Warning:can not make the fastq_2 file in this directory.\n");exit(1);}
		printf("the fastq_2 file:%s.\n",fastq);
	}

    fp=fopen(list,"r");							//open the species' list
	if(fp==NULL){printf("Warning:can not find the composition table file in this directory.\n");exit(1);}

    findex=fopen(index,"r");					//open the index file
	if(findex==NULL){printf("Warning:can not find the index file in this directory.\n");exit(1);}

    fc=fopen(config,"r");						//open the config file
	if(fc==NULL){printf("Warning:can not find the configure file in this directory.\n");exit(1);}

	sd_ratio[0]=2;

	//decide to simulate gemone or 16s;   0 for genome and 1 for 16s rRNA;
	int lable_type=0;

	if((strstr(simu_type,"16s")!=NULL)||(strstr(simu_type,"16S")!=NULL))
	{
	    lable_type=1;
	    printf("This is a 16s rRNA simulation!\n");
	    //totalcopy=count_total_copy(fp,findex,copy_number);
	    //printf("The total copy number is:%d\n",totalcopy);
	}

	else if(lable_type==0)
	{
	    printf("This is a genome simulation!\n");
	}
	else if(lable_type != 0)
	{
	    printf("This program can only simulate genome sequence or 16S rRNA sequence.\n");
	}
	//decide which method
	if(strstr(method,"illumina")!=NULL)
	{
		err_rand=(float *)malloc(200*50*sizeof(float));
		memset(err_rand,'\0',200*50*sizeof(float));
		min=30;  //reads' min length
		lable_method=2;
	}
	else if(strstr(method,"454")!=NULL)
	{
		err_rand=(float *)malloc(4000*50*sizeof(float));
		memset(err_rand,'\0',4000*50*sizeof(float));
		min=40;
		lable_method=1;
	}
	else if(strstr(method,"sanger")!=NULL)
	{
		length_max=length;
		min=0;
		lable_method=3;
	}
	else if(strstr(method,"pacbio")!=NULL)
	{
	    err_rand=(float *)malloc(50000*50*sizeof(float));
	    memset(err_rand,'\0',50000*50*sizeof(float));
	    min=150;
		lable_method=4;
		//length=450;
	}
	else
	{
		printf("This system only can simulate pacbio 454 illumina or sanger datas!\n");
		exit(1);
	}

	//decide error model;
	length_rand[0]=2; //as a flag
	if(lable_method!=3)
		length_max=get_qv_type(fc,type,subratio,sd_ratio,err_rand,lable_method,length_rand,exact);//function and get the simulation type
	input_right(length,reads,index,lable,gap_lable,gap,length_max,method,exact,length_rand);

	if(exact==0)
		if(sd_ratio[0]!=2)    //get the sd_ratio from config file
			sd=sd_ratio[0]*length;
		else if(lable_method==1)   //454,initiation
			sd=0.02*length;
		else
			sd=0;

	srand((unsigned)time(0));				//seed the rand
	total=species_number(fp,everynumber,name1,reads);				//function and get the whole number of species
	/*
	if(lable_type==1) //re_culculate the reads number according to the copy number
	{
		for(i=0;i<total;i++)
		{
		    printf("%s : %d \n",name1[i],everynumber[i]);
		}
		printf("After the reculculate by copy number:\n");
		re_everyread(everynumber,copy_number,total,reads);
	}
    */
	if(lable==1)						//malloc space for seed and lines
	{
		line2_1=(char *)malloc((length_max+1)*sizeof(char));
		line2_2=(char *)malloc((length_max+1)*sizeof(char));
		line4_1=(char *)malloc((length_max+1)*sizeof(char));
		line4_2=(char *)malloc((length_max+1)*sizeof(char));
	}
	else
	{
		line2=(char *)malloc((length_max+1)*sizeof(char));
		line4=(char *)malloc((length_max+1)*sizeof(char));
	}
	path=(char *)malloc(2000*sizeof(char));

	for(circle=0;circle<total;circle++)				 //in a specie,assign reads number
	{
        species_reads=everynumber[circle];
        printf("%s : %d\n",name1[circle],species_reads);

		flag=0;                         //flag:index has this species or not
		while(fgets(name,MAX_LINE,findex)!=NULL)					//get the specie's information in the index
		{
			name[strlen(name)-1]='\0';
    		split(index_name,name,"\t");
    		if((lable_type==0)&&(memcmp(name1[circle],index_name[3],strlen(name1[circle]))==0)&&(atoi(index_name[2])==0))  //add one condition
    		{
				flag=1;
				len=atoi(index_name[1]);
				memset(path,'\0',2000*sizeof(char));
				strcat(path,index_name[4]);
				gi=atoi(index_name[0]);
				break;
			}
			else if((lable_type==1)&&(memcmp(name1[circle],index_name[3],strlen(name1[circle]))==0))  //add one condition
    		{
				flag=1;
				len=atoi(index_name[1]);
				memset(path,'\0',2000*sizeof(char));
				strcat(path,index_name[4]);
				gi=atoi(index_name[0]);
				copy=atoi(index_name[2]);
				break;
			}
		}
		rewind(findex);

		if(flag==0)
		{
			printf("Warning:%s can't find in the index!!!",name1[circle]);
			continue;
		}

		//simulation
		buffle=(char *)malloc(len*sizeof(char));
		memset(buffle,'\0',len*sizeof(char));
		fpath=fopen(path,"r");
		if(fpath==NULL)
		{
			printf("Warning:can not find the index file pathway in this directory.\n");
			exit(1);
		}

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
		if(lable_type==1)
		{
		    minus_strand(buffle,minus_buffle,len);
		}
		else {minus_strand(buffle,minus_buffle,len);}


		flag_buff=0;           //buff's flag
		buff_point=0;		   //record the positon in the sequence
		buff_point1=0;
		buff_point2=0;
///////
		if(flag_bias==1)   //coverage bias
		{
			flag_gi=0;
			fbias=fopen(bias,"r");
			if(fbias==NULL){printf("Warning:can not find the coverage bias file.\n");exit(1);}
			all=(char *)malloc(NN*11*sizeof(char));
			memset(all,'\0',NN*11*sizeof(char));
			while(fgets(all,NN*11*sizeof(char),fbias)!=NULL)
			{
				all[strlen(all)-1]='\0';
				split(part_t,all,"=");

				get_gi=atoi(part_t[0]);
				if(gi==get_gi)
				{
					max_t=split(value_t,part_t[1],":");
					flag_gi=1;
					break;
				}
			}
			fclose(fbias);
			if(flag_gi==0)
			{
				printf("Warning:can not find the coverage bias information about %s\n",name1[circle]);
				exit(1);
			}
			for(x=0;x<max_t;x++)
				value_f[x]=atof(value_t[x]);
			free(all);
		}
///////////
		if(lable==1)			//simulation pair-end
		{
			char *buff_sequence1,*buff_sequence2;
			buff_sequence1=(char *)malloc((5*length_max+len_sequence_name)*buff);
			memset(buff_sequence1,'\0',(5*length_max+len_sequence_name)*buff);
			buff_sequence2=(char *)malloc((5*length_max+len_sequence_name)*buff);
			memset(buff_sequence2,'\0',(5*length_max+len_sequence_name)*buff);
			int every_gap;

			for(x=0;x<species_reads;x++)  //simulate everyread
			{
				memset(line2_1,'\0',length_max*sizeof(char));
				memset(line2_2,'\0',length_max*sizeof(char));
				memset(line4_1,'\0',length_max*sizeof(char));
				memset(line4_2,'\0',length_max*sizeof(char));
				if(exact==0) //use default length
				{
				    every_length=rand_length(length,sd,length_max,min);
				}
				else   //get length from input model
					every_length=length_distribution(length_rand,lable_method);
                every_gap=rand_length2(gap,'1',gap*10);
				if ((float)(rand())/(float)(RAND_MAX)<=0.5)                        //caculate the "+" "-" strain and the start position
				{
					flag_pn=0;
					if(flag_bias==1)
					{
						point1=start_position(max_t,flag_pn,len,(2*every_length+every_gap),value_f);
						point2=point1-2*every_length-every_gap;
					}
					else
					{
						point1=(int)(rand()/(float)(RAND_MAX)*(len-2*every_length-every_gap)+2*every_length+every_gap);
						point2=point1-2*every_length-every_gap;
					}
					buff_point1=first_line(point1,buff_sequence1,sequence_name,buff_point1,len_sequence_name,flag_pn);
					buff_point2=first_line(point2,buff_sequence2,sequence_name,buff_point2,len_sequence_name,(flag_pn*(-1)+1));
					buff_point1=bp(lable_type,lable_method,type,line2_1,line4_1,minus_buffle,subratio,buff_sequence1,buff_point1,every_length,err_rand,(len-point1-1),qual,len);
					buff_point2=bp(lable_type,lable_method,type,line2_2,line4_2,buffle,subratio,buff_sequence2,buff_point2,every_length,err_rand,point2,qual,len);
				}
				else
				{
					flag_pn=1;
					if(flag_bias==1)
					{
						point1=start_position(max_t,flag_pn,len,(2*every_length+every_gap),value_f);
						point2=point1+2*every_length+every_gap;
					}
					else
					{
						point1=(int)(rand()/(float)(RAND_MAX)*(len-2*every_length-every_gap));
						point2=point1+2*every_length+every_gap;
					}
					buff_point1=first_line(point1,buff_sequence1,sequence_name,buff_point1,len_sequence_name,flag_pn);
					buff_point2=first_line(point2,buff_sequence2,sequence_name,buff_point2,len_sequence_name,(flag_pn*(-1)+1));
					buff_point1=bp(lable_type,lable_method,type,line2_1,line4_1,buffle,subratio,buff_sequence1,buff_point1,every_length,err_rand,point1,qual,len);
					buff_point2=bp(lable_type,lable_method,type,line2_2,line4_2,minus_buffle,subratio,buff_sequence2,buff_point2,every_length,err_rand,(len-point2-1),qual,len);
				}

				buff_point1=remain_lines(buff_sequence1,buff_point1,line2_1,line4_1,every_length);     //function,add the ramain lines to sequence
				buff_point2=remain_lines(buff_sequence2,buff_point2,line2_2,line4_2,every_length);
				flag_buff++;

				if(flag_buff==buff||x==species_reads)     //write the sequence to hard disk
				{
					fprintf(fs1,"%s",buff_sequence1);
					fprintf(fs2,"%s",buff_sequence2);
					memset(buff_sequence1,'\0',(5*length_max+len_sequence_name)*buff);
					memset(buff_sequence2,'\0',(5*length_max+len_sequence_name)*buff);
					buff_point1=0;
					buff_point2=0;
					flag_buff=0;
				}
			}
			free(buff_sequence1);
			free(buff_sequence2);
		}
		else if(lable==0)											//simulation   single
		{
			char *buff_sequence;
			buff_sequence=(char *)malloc((5*length_max+len_sequence_name)*buff*sizeof(char));
			memset(buff_sequence,'\0',(5*length_max+len_sequence_name)*buff*sizeof(char));
			for(x=0;x<species_reads;x++)  //simulate everyread
			{
				memset(line2,'\0',length_max*sizeof(char));
				memset(line4,'\0',length_max*sizeof(char));
				if(exact==0)
					if(lable_method != 4)
				    {
				        every_length=rand_length(length,sd,length_max,min);
				    }
				    else
				    {
				        every_length=rand_length2(length,'0',length_max);
				    }
				else
					every_length=length_distribution(length_rand,lable_method);
				if ((float)(rand())/(float)(RAND_MAX)<=0.5)
				{
					flag_pn=0;
					if(flag_bias==1)
					{
						point=start_position(max_t,flag_pn,len,every_length,value_f);
						checkposition(flag_pn,point,len,every_length);
                    }
					else
						point=(int)(rand()/(float)(RAND_MAX)*(len-every_length)+every_length);
					buff_point=first_line(point,buff_sequence,sequence_name,buff_point,len_sequence_name,flag_pn);
					buff_point=bp(lable_type,lable_method,type,line2,line4,minus_buffle,subratio,buff_sequence,buff_point,every_length,err_rand,(len-point-1),qual,len);
				}
				else
				{
					flag_pn=1;
					if(flag_bias==1)
					{
						point=start_position(max_t,flag_pn,len,every_length,value_f);
						checkposition(flag_pn,point,len,every_length);
                    }
					else
						point=(int)(rand()/(float)(RAND_MAX)*(len-every_length));
					buff_point=first_line(point,buff_sequence,sequence_name,buff_point,len_sequence_name,flag_pn);
					buff_point=bp(lable_type,lable_method,type,line2,line4,buffle,subratio,buff_sequence,buff_point,every_length,err_rand,point,qual,len);
				}
				buff_point=remain_lines(buff_sequence,buff_point,line2,line4,every_length);
				flag_buff++;
				if(flag_buff==buff||x==species_reads)
				{
					fprintf(fs,"%s",buff_sequence);
					memset(buff_sequence,'\0',(5*length_max+len_sequence_name)*buff*sizeof(char));
					flag_buff=0;
					buff_point=0;
				}
			}
			free(buff_sequence);
		}
   		free(buffle);
		free(ref_name);
		free(sequence_name);
		free(minus_buffle);
	}
	fclose(fp);
	if(lable==1)
	{
		fclose(fs1);
		fclose(fs2);
		free(line2_1);
		free(line4_1);
		free(line2_2);
		free(line4_2);
	}
	else
	{
		fclose(fs);
		free(line2);
		free(line4);
	}
	fclose(findex);
	free(path);
	if(lable_method!=3)
		free(err_rand);
}

//======================================================================================
//function:buff_point_bp--add error position and type to result
//======================================================================================
int buff_point_bp(int position_err,int buff_point,char a1,char a2,char *buff_sequence)
{
	char s[6];
	sprintf(s,"%d",position_err);
	buff_sequence[buff_point]='|';
	buff_point++;
	strcat(&buff_sequence[buff_point],s);
	buff_point+=strlen(s);
	buff_sequence[buff_point]=':';
	buff_point++;
	buff_sequence[buff_point]=a1;
	buff_point++;
	buff_sequence[buff_point]=':';
	buff_point++;
	buff_sequence[buff_point]=a2;
	buff_point++;
	return buff_point;
}
//======================================================================================
//function:bp--simulate a read
//======================================================================================
int bp(int lable_type,int lable_method,float type[],char *line2,char *line4,char *buffle,float subratio[],char *buff_sequence,int buff_point,int length,float *err_rand,int position,float qual[],int len_seq)
{
	int a,turn,add_mark=0;
	float random,err_random,p;
	int quality();

	if(lable_method==3)  //sanger
	{
		for(a=0;a<length;a++)
		{
			line2[a]=buffle[a+position];
			line4[a]='r';
		}
	}
	else  //illumina,454,pacbio
	{
		for(a=0;(a+add_mark<length)&&((a+add_mark+position)<len_seq);a++)
		{
			random=(float)(rand())/(float)(RAND_MAX);
			turn=quality((a+add_mark),err_rand);
			p=qual[turn];
			if(random<=p)
   			{
				err_random=(float)(rand())/(float)(RAND_MAX);
				if(err_random<=type[0])  //substitute
				{

					if(lable_type==1){line2[a+add_mark]=sub_all(buffle[a+position],subratio[0],subratio[1],subratio[2],subratio[3],subratio[4],subratio[5],subratio[6],subratio[7],subratio[8],subratio[9],subratio[10],subratio[11]);}
					else{line2[a+add_mark]=sub_all(buffle[a+position],subratio[0],subratio[1],subratio[2],subratio[3],subratio[4],subratio[5],subratio[6],subratio[7],subratio[8],subratio[9],subratio[10],subratio[11]);}
					buff_point=buff_point_bp(a+add_mark,buff_point,buffle[a+position],line2[a+add_mark],buff_sequence);
					if(lable_method==2)
						line4[a+add_mark]=(char)(turn+33);
					else
						line4[a+add_mark]=(char)(turn+33);
				}
				else if(err_random<=type[1])   //insert
				{

					if(lable_type==1){line2[a+add_mark]=ins_all();}
					else{line2[a+add_mark]=ins_all();}
					buff_point=buff_point_bp(a+add_mark,buff_point,'-',line2[a+add_mark],buff_sequence);
					if(lable_method==2)
						line4[a+add_mark]=(char)(turn+33);
					else
						line4[a+add_mark]=(char)(turn+33);
					add_mark++;
					if((a+add_mark)<length)
					{
						line2[a+add_mark]=buffle[a+position];
						turn=quality((a+add_mark),err_rand);
						if(lable_method==2)
							line4[a+add_mark]=(char)(turn+33);
						else
							line4[a+add_mark]=(char)(turn+33);
					}
				}
				else    //delete
				{
					buff_point=buff_point_bp(a+add_mark,buff_point,buffle[a+position],'-',buff_sequence);
					add_mark--;
				}
			}
			else
			{
				line2[a+add_mark]=buffle[a+position];
				if(lable_method==2)
					line4[a+add_mark]=(char)(turn+33);
				else
					line4[a+add_mark]=(char)(turn+33);
			}
		}
	}
	return buff_point;
}
//======================================================================================
//function:input_right--judge the inputed parameters right or wrong
//======================================================================================
void input_right(int length,int reads,char index[],int lable,int gap_lable,int gap,int length_max,char method[],int exact,float length_rand[])
{
	printf("length=%d\n",length);
	if((length>length_max)&&(exact==0))
	{
		printf("Warning:The max read's length is:%d bps under the config file, so your reads' length can't exceed the max length!\n",length_max);
		exit(1);
	}

	printf("the reads number: %d.\n",reads);
	printf("the index file: %s.\n",index);
	printf("the method: %s.\n",method);
	if(exact==1)
	{
		if(length_rand[0]==2)
		{
			printf("this simulation-config-file can't allow exact simulation!\n");
			exit(1);
		}
		printf("this simulation is exact.\n");
	}
	else
		printf("the read length: %dbp.\n",length);
    printf("the max read length: %dbp.\n",length_max);
	if(lable==0)
	{
		if(gap_lable==1)
		{
			printf("Warning:You want to get the pair-end file,but this is the single model.Please use the parameter '-e 1'.\n");
			exit(1);
		}
		printf("this is the single model.\n");
	}
	if(lable==1)
	{
		printf("this is a pair-end model.\n");
		printf("the most width of pairs %d.\n",gap);
	}
}
//======================================================================================
//function:minus_strand--get the "-" strand according to the "+" strand
//======================================================================================
void minus_strand(char *buffle,char *minus_buffle,int len)
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
//======================================================================================
//function:minus_strand_r--get the "-" strand according to the "+" strand
//======================================================================================
void minus_strand_r(char *buffle,char *minus_buffle,int len)
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
//======================================================================================
//function:get_qv_type--get error config and return the simulation type
//======================================================================================
int get_qv_type(FILE *fc,float type[],float subratio[],float sd_ratio[],float *err_rand,int lable_method,float length_rand[],int exact)
{
	int circle_f,turn;
	char qv_config[40000];
	float value[72];
	int length_max=0;
	int get_value();
	turn=0;

	while(fgets(qv_config,40000,fc)!=NULL)
	{
		qv_config[strlen(qv_config)-1]='\0';
		if((strstr(qv_config,"illumina_rand")!=NULL) && (lable_method==2))  //get qv
		{
			get_value(qv_config,value);
			for(circle_f=0;circle_f<72;circle_f++)
				err_rand[length_max*72+circle_f]=value[circle_f];
			length_max++;
			continue;
		}
		else if((strstr(qv_config,"454_rand")!=NULL) && (lable_method==1))
		{
			get_value(qv_config,value);
			for(circle_f=0;circle_f<72;circle_f++)
				err_rand[length_max*72+circle_f]=value[circle_f];
			length_max++;
			continue;
		}
		else if((strstr(qv_config,"pacbio_rand")!=NULL) && (lable_method==4))
		{
			get_value(qv_config,value);
			for(circle_f=0;circle_f<72;circle_f++)
				err_rand[length_max*72+circle_f]=value[circle_f];
			length_max++;
			continue;
		}

		if((strstr(qv_config,"illumina_type")!=NULL && (lable_method==2)) ||(strstr(qv_config,"454_type_probability")!=NULL && (lable_method==1)) || (strstr(qv_config,"pacbio_type_probability")!=NULL && (lable_method==4)))
		{
			get_value(qv_config,type);
			continue;
		}

		if((strstr(qv_config,"454_sub")!=NULL && (lable_method==1))||(strstr(qv_config,"illumina_sub")!=NULL && (lable_method==2))||(strstr(qv_config,"pacbio_sub")!=NULL && (lable_method==4)))
		{
			get_value(qv_config,subratio);
			continue;
		}

		if((exact==0)&&(strstr(qv_config,"sd_ratio")!=NULL))
			get_value(qv_config,sd_ratio);
		else if((exact==1)&&(strstr(qv_config,"length_rand")!=NULL))
		{
			length_max=get_value(qv_config,length_rand);
			printf("The length_rand[0] is:%f",length_rand[0]);
			if(lable_method==1)
				length_max=length_max+39;
			else if(lable_method==2)
				length_max=length_max+29;
            else if(lable_method==4)
				length_max=length_max+149;
		}
	}
	fclose(fc);
	return length_max;
}
//======================================================================================
//function:get_value
//======================================================================================
int get_value(char qv_config[],float need[])
{
	int circle_f,length;
	char *part[2],*value[40000];
	split(part,qv_config,"=");
	length=split(value,part[1],":");
	for(circle_f=0;circle_f<length;circle_f++)
	{
		need[circle_f]=atof(value[circle_f]);
	}
	//printf("%f%f%f\n",need[0],need[1],need[2]);
	return length;
}
//======================================================================================
//function:species_number--how many reads in every species and return how many species
//======================================================================================
int species_number(FILE *fp,int everynumber[],char name1[2000][200],int reads)
{
	int total=0;
	float percentage1[2000][1];
	float phase[2000],random;
	int circle1,circle2;
	char species_percent[MAX_LINE],*species_name[MAX_LINE];
	while(fgets(species_percent,MAX_LINE,fp)!=NULL)
	{
		species_percent[strlen(species_percent)-1]='\0';
		split(species_name,species_percent,"\t");
		strcpy(name1[total],species_name[0]);
		percentage1[total][0]=atof(species_name[1]);
		everynumber[total]=0;
		total++;
	}
	for(circle1=0;circle1<total-1;circle1++)
	{
		phase[circle1]=0;
		for(circle2=0;circle2<=circle1;circle2++)
		{
			phase[circle1]=phase[circle1]+percentage1[circle2][0];
		}
	}
	for(circle1=0;circle1<reads;circle1++)
	{
		random=(float)rand()/(float)(RAND_MAX);
		if(random>phase[total-2])
		{
			everynumber[total-1]++;
		}
		else if(random<=phase[0])
		{
			everynumber[0]++;
		}
		else
		{
			for(circle2=0;circle2<total-2;circle2++)
			{
				if(random>phase[circle2]&&random<=phase[circle2+1])
				{
					everynumber[circle2+1]++;
					continue;
				}
			}
		}
	}
	return total;
}
//======================================================================================
//function:first_line--add the reference name,start position and "+" "-"strain to the sequence
//======================================================================================
int first_line(int point,char *buff_sequence,char *sequence_name,int buff_point,int len_sequence_name,int flag_pn)//flag=1:the first in pair-end and 2
{
	char prints[20];
	strcat(&buff_sequence[buff_point],sequence_name);
	buff_point=buff_point+len_sequence_name;
	buff_sequence[buff_point]='|';
	buff_point++;
	sprintf(prints,"%d",point);
	strcat(&buff_sequence[buff_point],prints);
	buff_point=buff_point+strlen(prints);
	buff_sequence[buff_point]='|';
	if (flag_pn==1)
		buff_sequence[buff_point+1]='+';
	else
		buff_sequence[buff_point+1]='-';
	buff_point+=2;
	return buff_point;
}
//======================================================================================
//function:remain_lines--add simulated sequence and it QV to the result:buff_sequence
//======================================================================================
int remain_lines(char *buff_sequence,int buff_point,char *line2,char *line4,int every_length)
{
	buff_sequence[buff_point]='\n';
	buff_point++;
	strcat(&buff_sequence[buff_point],line2);
	buff_point+=every_length;
	buff_sequence[buff_point]='\n';
	buff_sequence[buff_point+1]='+';
	buff_sequence[buff_point+2]='\n';
	buff_point+=3;
	strcat(&buff_sequence[buff_point],line4);
	buff_point+=every_length;
	buff_sequence[buff_point]='\n';
	buff_point++;
	return buff_point;
}
//======================================================================================
//function:quality--return an error value
//======================================================================================
int quality(int pos_base,float *err_rand)  //pos_base:from zero
{
	float random;
	int circle,step,turn;
	step=72;
	random=(float)(rand())/(float)(RAND_MAX);

	if(random<=err_rand[step*pos_base+30])
	{
		for(circle=(step*pos_base+30);circle>step*pos_base;circle--)
		{
			if((random>err_rand[circle-1])&&(random<=err_rand[circle]))
			{
				turn=circle-step*pos_base;
				return turn;
			}
		}
		return 0;
	}
	else
		for(circle=(step*pos_base+31);circle<step*(pos_base+1);circle++)
		{
			if((random>err_rand[circle-1])&&(random<=err_rand[circle]))
			{
				turn=circle-step*pos_base;
				return turn;
			}
		}
}
//======================================================================================
//function:length_distribution--return an random length
//======================================================================================
int length_distribution(float length_rand[],int lable_method)
{
	float random;
	int i,length;
	random=(float)(rand())/(float)(RAND_MAX);
	if(random<=length_rand[0])
	{
		if(lable_method==1)
			length=40;
		else if(lable_method==2)
			length=30;
        else if(lable_method==4)
            length=150;
		return length;
	}
	else
	{
		for(i=0; ;i++)
		{
			if((random>length_rand[i])&&(random<=length_rand[i+1]))
			{
				if(lable_method==1)
					length=41+i;
				else if(lable_method==2)
					length=31+i;
                else if(lable_method==4)
                    length=150+i;
				return length;
			}
		}
	}
}
//======================================================================================
//function:start_position--return an random start position
//======================================================================================
int start_position(int max_t,int flag_pn,int len,int every_length,float value_f[])
{
	int i,pos,flag,start;
	float random;

	flag=0;
	pos=-1;
	while(flag==0)
	{
		random=(float)(rand())/(float)(RAND_MAX);
		start=(int)(max_t*random);
		if(random>value_f[start])
		{
			for(i=start;i<max_t;i++)
			{
				if(random<value_f[i])
				{
					pos=i;
					break;
				}
			}
		}
		else
		{
			for(i=start;i>=0;i--)
			{
				if(random>value_f[i])
				{
					pos=i+1;
					break;
				}
			}
			if(pos==-1)
			{
				pos=0;
				break;
			}
		}

		pos=pos*100+rand()%100;
		if(flag_pn==1)
		{
			if((len-every_length)>=pos)
				flag=1;
			else
				flag=0;
		}
		else
		{
			if((pos>=every_length)&&(pos<len))
				flag=1;
			else
				flag=0;
		}
	}
	return pos;
}

void checkposition(int flag_pn,int point,int len, int everylength)
{
    if(flag_pn==0)
    {
        if((point + everylength) > len)
        printf("ERROR READ (forward)!!!\n%d\t%d\t%d\n",point,len,everylength);
    }
    else
    {
        if((point < everylength)||(point >= len))
        printf("ERROR READ (revised)!!!\n%d\t%d\t%d\n",point,len,everylength);
    }

}
//======================================================================================
//function:count_total_copy--return the total copy number of the specices in the composition table
//======================================================================================
int count_total_copy(FILE *fp,FILE *findex,int copy_number[])
{
    int total_copy=0;
    int copy=0,total=0;
	int circle;
	char species_percent[MAX_LINE],*species_name[MAX_LINE];
	char name1[2000][200];
	char name[MAX_LINE],*index_name[MAX_LINE];
    while(fgets(species_percent,MAX_LINE,fp)!=NULL)
    {
        species_percent[strlen(species_percent)-1]='\0';
		split(species_name,species_percent,"\t");
		strcpy(name1[total],species_name[0]);
		total++;
    }
    rewind(fp);
    for(circle=0;circle<total;circle++)
    {
        while(fgets(name,MAX_LINE,findex)!=NULL)					//get the specie's information in the index
		{
			name[strlen(name)-1]='\0';
    		split(index_name,name,"\t");
			if((memcmp(name1[circle],index_name[3],strlen(name1[circle]))==0))  //add one condition
    		{
				copy_number[circle]=atoi(index_name[2]);
				printf("%s copynumber is %d\n",name1[circle],copy_number[circle]);
				total_copy=copy_number[circle]+total_copy;
				break;
			}
		}
		rewind(findex);
    }
    return total_copy;
}
//======================================================================================
//function:re_everyread--return the everyread number after considering the copynumber.
//======================================================================================
void re_everyread(int everyread[],int copy_number[],int total,int reads)
{
    int tmp[2000];
    int new_total=0;
    int i;
    for(i=0;i<total;i++)
    {

        tmp[i]=everyread[i]*copy_number[i];

        new_total=new_total+tmp[i];
    }
    for(i=0;i<total;i++)
    {
        everyread[i]=reads*tmp[i]/new_total;

    }
}
