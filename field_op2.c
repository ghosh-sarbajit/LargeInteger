#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

/*=============================================================================
Objective:This program is for calculation of modular reduction of over the prime
field p(=2^127-1), also using modular reduction we will calculate finite field
arithmetic, e.g. field addition, multiplication, subtraction, inverse.

Input: Two integers, in the range {0, 1, ..., 2^127-2} in hexadecimal.

Output: i. Addition, subtraction, multiplication of two field elements.
	ii. Also  inverse of both the elements in the field.

Modular reduction is implemented with in both addition and multiplication func.

We have implemented the the inverse function without using any if, else condition.
=============================================================================*/

/* To print dashed line for formatted output*/
void line_print(int k)
{
    int j;
    for(j=0;j<k;j++) putchar('*');
    putchar('\n');
}

//polynomial representation of two integers using 5 box of 64 bit and linmb size is 28 bit
uint64_t poly_1[5]={0UL, 0UL, 0UL, 0UL, 0UL};
uint64_t poly_2[5]={0UL, 0UL, 0UL, 0UL, 0UL};
//polynmial representation of two integes after sum and mult
uint64_t sum[5]={0UL, 0UL, 0UL, 0UL, 0UL};
uint64_t diff[5]={0UL, 0UL, 0UL, 0UL, 0UL};
uint64_t mult[9]={0UL, 0UL, 0UL, 0UL, 0UL, 0UL, 0UL, 0UL, 0UL}; //mult of two 4 deg poly will yield a 8 deg poly
uint64_t mult_high[5]={0UL, 0UL, 0UL, 0UL, 0UL}; //upper 127 bit of mult poly
uint64_t mult_low[5]={0UL, 0UL, 0UL, 0UL, 0UL}; //lower 127-bit of mult poly
uint64_t mult_poly[5]={0UL, 0UL, 0UL, 0UL, 0UL}; // used to store poly representaion of field multiplication
uint64_t inv_result_poly[5]={0UL, 0UL, 0UL, 0UL, 0UL};
//return from poly representation
uint64_t sum_result[2]={0UL,0UL};
uint64_t num_2_two_c[2]={0UL,0UL};
uint64_t diff_result[2]={0UL,0UL};
uint64_t mult_result[2]={0UL,0UL};
uint64_t inv_result[2]={0UL, 0UL};
//calculation of 2's complement
uint64_t poly_2_one_c[5]={0UL, 0UL, 0UL, 0UL, 0UL}; /*1's complement of poly_2 will be saved here*/
uint64_t poly_2_two_c[5]={0UL, 0UL, 0UL, 0UL, 0UL}; /*2's complement of poly_2 will be saved here*/
uint64_t one[5]={0x1, 0UL, 0UL, 0UL, 0UL};
uint64_t _p[5]={0xfffffff, 0xfffffff, 0xfffffff, 0xfffffff, 0x7fff};

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

/*==============================================================================
addition function takes 4 arguments
first 2 args are pointer to array of poly represention of two number
3rd arg is pointer to the array where it stores poly representation of their sum
4th arg is ptr to the array where it stores int representation of their sum, using two two 64 bit int
===============================================================================*/
void field_add_two_num(uint64_t *poly_1, uint64_t *poly_2, uint64_t *sum, uint64_t *sum_result)
{
	uint64_t carry=0, temp_num=0;
    for(int count1=0; count1<5; count1++)   sum[count1]=0UL;
    for(int count1=0; count1<2; count1++)   sum_result[count1]=0UL;
	for(int count1=0; count1<5;count1++)
	{
		temp_num=0;
		temp_num= poly_1[count1]+poly_2[count1]+carry;
		sum[count1]=temp_num&0xfffffff;
		carry=temp_num>>28;
	}
    /*_128_bit will be used to check whether overflow has ocurred*/
    uint64_t _128_bit=sum[4]>>15;
    sum[4]=sum[4]&0x7fff;
    /* Adjusting sum according to overflow*/
    /* observation: if we add (p-1)+(p-1) then last two bytes are 0*/
    /* So we can always accomodate 1-bit at last pos, or the is of 127 bit*/
    temp_num=(sum[0]+_128_bit)*(_128_bit==1)+(sum[0])*(_128_bit!=1);
    sum[0]=temp_num&0xfffffff;

    /* Note p(=pow(2,127)) is 0 in F_p*/
    uint64_t indicator;
    indicator=(sum[0]==0xfffffff)&&(sum[1]==0xfffffff)&&(sum[2]==0xfffffff)&&(sum[3]==0xfffffff)&&(sum[4]==0x7fff);
    sum[0]=(sum[0])*(indicator!=1)+(0x0)*(indicator==1);
    sum[1]=(sum[1])*(indicator!=1)+(0x0)*(indicator==1);
    sum[2]=(sum[2])*(indicator!=1)+(0x0)*(indicator==1);
    sum[3]=(sum[3])*(indicator!=1)+(0x0)*(indicator==1);
    sum[4]=(sum[4])*(indicator!=1)+(0x0)*(indicator==1);

	sum_result[0]=sum_result[0]|sum[0]; //7
	sum_result[0]=sum_result[0]|(sum[1]<<28); //7+7
	sum_result[0]=sum_result[0]|(sum[2]<<56); // 7+7+2 = 16
	sum_result[1]=sum_result[1]|(sum[2]>>8); // 5 // lower 2 byte consumed earlier
	sum_result[1]=sum_result[1]|(sum[3]<<20); // 5+7
	sum_result[1]=sum_result[1]|(sum[4]<<48); // 5+7+4
}


/*=============================================================================
multiplicatin function
multiplication takes 4 arguments
first two are ptr to the array of the poly representation of two integers
3rd is the ptr where it saves the value of field mult in poly format
4th is the ptr where it saves the value of field mult in compact form
=============================================================================*/
void field_mult_two_num(uint64_t *poly_1, uint64_t *poly_2, uint64_t *mult_poly, uint64_t *mult_result)
{
	uint64_t carry=0, temp_num=0;

	temp_num=poly_1[0]*poly_2[0];
	mult[0]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[0]*poly_2[1]+poly_1[1]*poly_2[0];
	temp_num=temp_num+carry;
	mult[1]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[0]*poly_2[2]+poly_1[1]*poly_2[1]+poly_1[2]*poly_2[0];
	temp_num=temp_num+carry;
	mult[2]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[0]*poly_2[3]+poly_1[1]*poly_2[2]+poly_1[2]*poly_2[1]+poly_1[3]*poly_2[0];
	temp_num=temp_num+carry;
	mult[3]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[0]*poly_2[4]+poly_1[1]*poly_2[3]+poly_1[2]*poly_2[2]+poly_1[3]*poly_2[1]+poly_1[4]*poly_2[0];
	temp_num=temp_num+carry;
	mult[4]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[1]*poly_2[4]+poly_1[2]*poly_2[3]+poly_1[3]*poly_2[2]+poly_1[4]*poly_2[1];
	temp_num=temp_num+carry;
	mult[5]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[2]*poly_2[4]+poly_1[3]*poly_2[3]+poly_1[4]*poly_2[2];
	temp_num=temp_num+carry;
	mult[6]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[3]*poly_2[4]+poly_1[4]*poly_2[3];
	temp_num=temp_num+carry;
	mult[7]=temp_num&0xfffffff;
	carry=temp_num>>28;

	temp_num=poly_1[4]*poly_2[4];
	temp_num=temp_num+carry;
	mult[8]=temp_num;

    /* making separation of lowe 127-bit*/
    mult_low[0]=mult[0];
    mult_low[1]=mult[1];
    mult_low[2]=mult[2];
    mult_low[3]=mult[3];
    mult_low[4]=mult[4]&0x7fff;

    mult_high[0]=mult[4]>>15;//13
    mult_high[0]=mult_high[0]|(mult[5]<<13);//13+15
    mult_high[0]&=0xfffffff;

    mult_high[1]=mult[5]>>15;//13
    mult_high[1]=mult_high[1]|(mult[6]<<13);//13+15
    mult_high[1]&=0xfffffff;

    mult_high[2]=mult[6]>>15;//13
    mult_high[2]=mult_high[2]|(mult[7]<<13);//13+15
    mult_high[2]&=0xfffffff;

    mult_high[3]=mult[7]>>15;//13
    mult_high[3]=mult_high[3]|(mult[8]<<13);//13+15
    mult_high[3]&=0xfffffff;

    mult_high[4]=mult[8]>>15;
    mult_high[4]&=0x7fff;

    field_add_two_num(mult_low, mult_high, mult_poly, mult_result);
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
1st arg is pointer to array of poly-represention of 1st num
2nd arg is pointer to array of 2's comp poly-represention of 2nd num
3rd arg is pointer to the array where it stores poly representation of their difference
4th arg is pointer to the array where it stores integer representation of their sum
==============================================================================*/
void diff_two_num(uint64_t *poly_1, uint64_t *poly_2, uint64_t *diff, uint64_t *diff_result)
{
	uint64_t carry=0, temp_num=0, indicator=0;// since we already passing complement
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
		find_twos_complement(diff, diff, diff, diff_result);
        for(int count1=0; count1<5; count1++) poly_2[count1]=diff[count1];
        find_twos_complement(poly_2, poly_2_one_c, poly_2_two_c, diff_result);
        diff_result[0]=0, diff_result[1]=0;
        diff_two_num(_p, poly_2_two_c, diff, diff_result);
        return;
	}
}
/*=============================================================================
Note over F_p, Inv(a)=pow(a, p-2), and p=pow(2,127)-1
Here we will use square and multiply method for this
Since here p is fixed, we have used the binary representation of p-2, which is
'0b1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111101'
so clearly we have to sq & multiply 125 times, and then one square and one square and multiply
===============================================================================*/
void field_inv(uint64_t *poly_1, uint64_t *inv_result_poly, uint64_t *inv_result)
{
    uint64_t temp_poly[5]={0UL, 0UL, 0UL, 0UL, 0UL};
    uint64_t temp_inv_result[2]={0UL, 0UL};
    // first iteration done (1^2)*a= a
    for(int count1=0; count1<5; count1++)   inv_result_poly[count1]=poly_1[count1];

    // from 2 to 125 th iteration
    for(int count2=2; count2<=125;count2++)
    {
        field_mult_two_num(inv_result_poly, inv_result_poly, inv_result_poly, inv_result);
        for(int count1=0; count1<5; count1++)   temp_poly[count1]=inv_result_poly[count1];
        for(int count1=0; count1<2; count1++)   temp_inv_result[count1]=inv_result[count1];
        field_mult_two_num(temp_poly, poly_1, inv_result_poly, inv_result);
        for(int count1=0; count1<5; count1++)   temp_poly[count1]=inv_result_poly[count1];
        for(int count1=0; count1<2; count1++)   temp_inv_result[count1]=inv_result[count1];
    }
    //126 th bit is zero, so only Square
    field_mult_two_num(inv_result_poly, inv_result_poly, inv_result_poly, inv_result);
    for(int count1=0; count1<5; count1++)   temp_poly[count1]=inv_result_poly[count1];
    for(int count1=0; count1<2; count1++)   temp_inv_result[count1]=inv_result[count1];
    //127 th bit is 1, so we will do square and multiply both
    field_mult_two_num(inv_result_poly, inv_result_poly, inv_result_poly, inv_result);
    for(int count1=0; count1<5; count1++)   temp_poly[count1]=inv_result_poly[count1];
    for(int count1=0; count1<2; count1++)   temp_inv_result[count1]=inv_result[count1];
    field_mult_two_num(temp_poly, poly_1, inv_result_poly, inv_result);
    for(int count1=0; count1<5; count1++)   temp_poly[count1]=inv_result_poly[count1];
    for(int count1=0; count1<2; count1++)   temp_inv_result[count1]=inv_result[count1];
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
	uint64_t temp_num=0, indicator;
	unsigned char ch;


	// taking user input for first integer
	line_print(60);
    printf("%s", "Program for calculation of finite field arithmatic\n");
    printf("%s", "Over the prime field of size p=pow(2,127)-1\n");
    printf("%s", "Please enter two field element in the range \n{0x0, 0x1, ..., p-1=0x7ffffffffffffffffffffffffffffffe}\n");
    line_print(60);
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

	field_add_two_num(poly_1, poly_2, sum, sum_result);
	printf("Field addition of two elements over F_p where p=pow(2,127)-1\n");
	for(int count1=1;count1>=0;count1--)
	{
		printf("%.16lx\t", sum_result[count1]);
	}
	printf("\n\n" );

	field_mult_two_num(poly_1, poly_2, mult_poly, mult_result);
	printf("Field multiplication of two elements over F_p where p=pow(2,127)-1\n");
	for(int count1=1;count1>=0;count1--)
	{
		printf("%.16lx\t", mult_result[count1]);
	}
	printf("\n\n" );

    //invese of entered elements
    field_inv(poly_1, inv_result_poly, inv_result);
    printf("Inverse of first element\n");
	for(int count1=1;count1>=0;count1--)
	{
		printf("%.16lx\t", inv_result[count1]);
	}
	printf("\n\n" );
    field_inv(poly_2, inv_result_poly, inv_result);
    printf("Inverse of second element\n");
	for(int count1=1;count1>=0;count1--)
	{
		printf("%.16lx\t", inv_result[count1]);
	}
	printf("\n\n" );

	// substration of two values
	find_twos_complement(poly_2, poly_2_one_c, poly_2_two_c, num_2_two_c);
	diff_two_num(poly_1, poly_2_two_c, diff, diff_result);


	printf("Field minus of two elements over F_p where p=pow(2,127)-1\n");
    printf("Note we have implemented as -a=p-a (mod p)\n" );
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
python3 -c 'print(hex( 0x7ffffffffffffffffffffffffffffffa * 0x7fffffffffffffffffffffffffffff ))'
3fffffffffffffffffffffffffffffff00000000000000000000000000000001
python3 -c 'print(hex( 0x2468b59754249e - 0x2468b59754 ))'
==============================================================================*/
