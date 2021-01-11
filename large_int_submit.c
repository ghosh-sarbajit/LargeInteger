#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

/*=============================================================================
This program is for calculation of 127-bit integer addition, substration, and
multiplicatin. It if two 127 bit integer is inserted for addition it can handle,
one bit overflow, in case of addition.
=============================================================================*/

//polynomial representation of two integers using 5 box of 64 bit and linmb size is 28 bit
uint64_t poly_1[5]={0UL, 0UL, 0UL, 0UL, 0UL};
uint64_t poly_2[5]={0UL, 0UL, 0UL, 0UL, 0UL};
//polynmial representation of two integes after sum and mult
uint64_t sum[5]={0UL, 0UL, 0UL, 0UL, 0UL};
uint64_t diff[5]={0UL, 0UL, 0UL, 0UL, 0UL};
uint64_t mult[9]={0UL, 0UL, 0UL, 0UL, 0UL, 0UL, 0UL, 0UL, 0UL}; //mult of two 4 deg poly will yield a 9 deg poly
//return from poly representation
uint64_t sum_result[2]={0UL,0UL};
uint64_t num_2_two_c[2]={0UL,0UL};
uint64_t diff_result[2]={0UL,0UL};
uint64_t mult_result[4]={0UL,0UL,0UL,0UL};
//calculation of 2's complement
uint64_t poly_2_one_c[5]={0UL, 0UL, 0UL, 0UL, 0UL}; /*1's complement of poly_2 will be saved here*/
uint64_t poly_2_two_c[5]={0UL, 0UL, 0UL, 0UL, 0UL}; /*2's complement of poly_2 will be saved here*/
uint64_t one[5]={0x1, 0UL, 0UL, 0UL, 0UL};

/*==============================================================================
addition function takes 4 arguments
first 2 args are pointer to array of poly represention of two number
3rd arg is pointer to the array where it stores poly representation of their sum
4th arg is ptr to the array where it stores int representation of their sum, using two two 64 bit int
===============================================================================*/
void add_two_num(uint64_t *poly_1, uint64_t *poly_2, uint64_t *sum, uint64_t *sum_result)
{
	uint64_t carry=0, temp_num=0;
	for(int count1=0; count1<5;count1++)
	{
		temp_num=0;
		temp_num= poly_1[count1]+poly_2[count1]+carry;
		sum[count1]=temp_num&0xfffffff;
		carry=temp_num>>28;
	}
	sum_result[0]=sum_result[0]|sum[0]; //7
	sum_result[0]=sum_result[0]|(sum[1]<<28); //7+7
	sum_result[0]=sum_result[0]|(sum[2]<<56); // 7+7+2 = 16
	sum_result[1]=sum_result[1]|(sum[2]>>8); // 5 // lower 2 byte consumed earlier
	sum_result[1]=sum_result[1]|(sum[3]<<20); // 5+7
	sum_result[1]=sum_result[1]|(sum[4]<<48); // 5+7+4
}

/*=============================================================================
multiplicatin function
multiplication takes two arguments
these are pointer array of the poly representation of two integers
=============================================================================*/
void mult_two_num(uint64_t *poly_1, uint64_t *poly_2)
{
	//multiplication of two integers
	// m_i = sum(p_j q_k) s.t i=j+k && 0<=j<5 && 0<=k<5, for all 0<=i<9
	// m_i = sum( p_j q_(i-j) ) s.t 0<=j<5, for all 0<=i<9

	uint64_t carry=0, temp_num=0;

	temp_num=poly_1[0]*poly_2[0];
	mult[0]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[0]*poly_2[1]+poly_1[1]*poly_2[0];
	temp_num=temp_num+carry;
	mult[1]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[0]*poly_2[2]+poly_1[1]*poly_2[1]+poly_1[2]*poly_2[0];
	temp_num=temp_num+carry;
	mult[2]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[0]*poly_2[3]+poly_1[1]*poly_2[2]+poly_1[2]*poly_2[1]+poly_1[3]*poly_2[0];
	temp_num=temp_num+carry;
	mult[3]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[0]*poly_2[4]+poly_1[1]*poly_2[3]+poly_1[2]*poly_2[2]+poly_1[3]*poly_2[1]+poly_1[4]*poly_2[0];
	temp_num=temp_num+carry;
	mult[4]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[1]*poly_2[4]+poly_1[2]*poly_2[3]+poly_1[3]*poly_2[2]+poly_1[4]*poly_2[1];
	temp_num=temp_num+carry;
	mult[5]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[2]*poly_2[4]+poly_1[3]*poly_2[3]+poly_1[4]*poly_2[2];
	temp_num=temp_num+carry;
	mult[6]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[3]*poly_2[4]+poly_1[4]*poly_2[3];
	temp_num=temp_num+carry;
	mult[7]=temp_num&0xfffffff;
	carry=temp_num>>28;
	temp_num=0;

	temp_num=poly_1[4]*poly_2[4];
	temp_num=temp_num+carry;
	mult[8]=temp_num;

	mult_result[0]=mult_result[0]|mult[0];//7
	mult_result[0]=mult_result[0]|(mult[1]<<28); //7+7
	mult_result[0]=mult_result[0]|(mult[2]<<56); //7+7+2
	mult_result[1]=mult_result[1]|(mult[2]>>8); //5
	mult_result[1]=mult_result[1]|(mult[3]<<20); // 5+7
	mult_result[1]=mult_result[1]|(mult[4]<<48); //5+7+4
	mult_result[2]=mult_result[2]|(mult[4]>>16);//3
	mult_result[2]=mult_result[2]|(mult[5]<<12);//3+7
	mult_result[2]=mult_result[2]|(mult[6]<<40);//3+7+6
	mult_result[3]=mult_result[3]|(mult[6]>>24);//1
	mult_result[3]=mult_result[3]|(mult[7]<<4);//1+7
	mult_result[3]=mult_result[3]|(mult[8]<<28);//1+7+7
	mult_result[3]=mult_result[3]|(mult[9]<<60);//1+7+7+1
}

/*=============================================================================
calculation of 2's complement
it takes 4 args
first is the poly representation of integer whose 2's comp we calculate
2nd is the pointer to array of integer where it stores 1's complement
3rd is the pointer to the array where it stores 2's complement, in poly representation
4th is the pointer to the array whrer it sores 2's complement, in integer representaion
==============================================================================*/
void find_twos_complement(uint64_t *source, uint64_t *intermideate, uint64_t *target1, uint64_t *target2)
{
	uint64_t temp_num;
	for(int count1=0; count1<4;count1++)
	{
		temp_num=0;
		temp_num=source[count1]^0xfffffff;
		intermideate[count1]=temp_num;
	}
	intermideate[4]=source[4]^0x7fff;
	add_two_num(intermideate, one, target1, target2);
}

/*=============================================================================
substraction function takes 4 arguments
1st arg is pointer to array of poly represention of 1st num
2nd arg is pointer to array of 2's com poly represention of 2nd num
3rd arg is pointer to the array where it stores poly representation of their difference
4th arg is pointer to the array where it stores integer representation of their sum
==============================================================================*/
void diff_two_num(uint64_t *poly_1, uint64_t *poly_2, uint64_t *diff, uint64_t *diff_result)
{
	uint64_t carry=0, temp_num=0, indicator=0;
	for(int count1=0; count1<5;count1++)
	{
		temp_num=0;
		temp_num= poly_1[count1]+poly_2[count1]+carry;
		diff[count1]=temp_num&0xfffffff;
		carry=temp_num>>28;
	}
	indicator=(diff[4]>>15);
	if(indicator==1)
	{
		diff[4]=diff[4]&0x7fff;
		diff_result[0]=diff_result[0]|diff[0]; //7
		diff_result[0]=diff_result[0]|(diff[1]<<28); //7+7
		diff_result[0]=diff_result[0]|(diff[2]<<56); // 7+7+2 = 16
		diff_result[1]=diff_result[1]|(diff[2]>>8); // 5 // lower 2 byte consumed earlier
		diff_result[1]=diff_result[1]|(diff[3]<<20); // 5+7
		diff_result[1]=diff_result[1]|(diff[4]<<48); // 5+7+4
		return;
	}
	else
	{
		printf("2nd integer is larger\n" );
		find_twos_complement(diff, diff, diff, diff_result);
		return;
	}
}


int main()
{
	//for taking user input corresponding to num1 and num2
	unsigned char *num_1, *num_2;
	num_1=malloc(32*sizeof(unsigned char));
	num_2=malloc(32*sizeof(unsigned char));
	for(int count1=0; count1<32;count1++) *(num_1+count1)=0;
	for(int count1=0; count1<32;count1++) *(num_2+count1)=0;

	// some intermideate temp vars
	uint64_t temp_num=0;
	unsigned char ch;


	// taking user input for first integer
	printf("%s", "Enter number one:");
	unsigned char i=0; //for tracking input length
	while((ch=getchar())!='\n')
	{
		if(ch >= 'a' && ch <= 'f')
		{
			num_1[i]=(ch-'a')+10;
		}
		if(ch >= 'A' && ch <= 'F')
		{
			num_1[i]=(ch-'A')+10;
		}
		if(ch >= '0' && ch <= '9')
		{
			num_1[i]=(ch-'0');
		}
		i++;
	}
	i--;

	//creating polynomial representation of first integer
	for(int count1=0; count1<5; count1++)
	{
		for(int count2=0; count2<7 && i>=0; count2++)
		{
			temp_num=(uint64_t)num_1[i];
			poly_1[count1]=poly_1[count1]|(temp_num<<(4*count2));
			i--;
		}
	}


	// taking user input for 2nd integer
	printf("%s", "Enter number two:");
	i=0;
	while((ch=getchar())!='\n')
	{
		if(ch >= 'a' && ch <= 'f')
		{
			num_2[i]=(ch-'a')+10;
		}
		if(ch >= 'A' && ch <= 'F')
		{
			num_2[i]=(ch-'A')+10;
		}
		if(ch >= '0' && ch <= '9')
		{
			num_2[i]=(ch-'0');
		}
		i++;
	}
	i--;

	//creating polynomial representation of second integer
	for(int count1=0; count1<5; count1++)
	{
		for(int count2=0; count2<7 && i>=0; count2++)
		{
			temp_num=(uint64_t)num_2[i];
			poly_2[count1]=poly_2[count1]|(temp_num<<(4*count2));
			i--;
		}
	}


	add_two_num(poly_1, poly_2, sum, sum_result);
	printf("Compact representation of sum\n");
	for(int count1=1;count1>=0;count1--)
	{
		printf("%.16lx\t", sum_result[count1]);
	}
	printf("\n" );
	mult_two_num(poly_1,poly_2);
	printf("Compact representation of mult\n");
	for(int count1=3;count1>=0;count1--)
	{
		printf("%.16lx\t", mult_result[count1]);
	}
	printf("\n" );

	// substration of two values
	find_twos_complement(poly_2, poly_2_one_c, poly_2_two_c, num_2_two_c);
	diff_two_num(poly_1, poly_2_two_c, diff, diff_result);

	printf("Compact representation of diff\n");
	for(int count1=1;count1>=0;count1--)
	{
		printf("%.16lx\t", diff_result[count1]);
	}
	printf("\n" );

	return 0;
}

/*=============================================================================
Sample run
0x2468b59754249e
0x1234aabb5678ccdd1234ffff1234eeee
python3 -c 'print(hex( 0x89a2359acdbe977 * 0x458907236889 ))'
python3 -c 'print(hex( 0x89a2359acdb123abdee977 * 0x4589073677236889 ))'
12345acbaa124f
python3 -c 'print(hex( 0x7fffffffffffffffffffffffffffffff * 0x4589073677236889 ))'
python3 -c 'print(hex( 0x7fffffffffffffffffffffffffffffff * 0x7fffffffffffffffffffffffffffffff ))'
3fffffffffffffffffffffffffffffff00000000000000000000000000000001
python3 -c 'print(hex( 0x2468b59754249e - 0x2468b59754 ))'
==============================================================================*/
