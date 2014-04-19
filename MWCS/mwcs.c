/* ******************
 * ** Assignment 3 **
 * ******************
 * 
 * File:   mwcs.c
 * Author: Pijus Kumar Sarker
 * Student ID: 103945877
 * 
 * Course: Algorithms in Bioinformatics, Fall 2013
 * Submitted: November 28, 2013
 *
 * Created on November 23, 2013, 12:53 PM
 * 
 * +---------------------------------------------------------------------------------------------------+
 * | PROBLEM                                                                                           |
 * | --------                                                                                          |
 * | The Mutation Sensitive Alignment (MSA) algorithm computes in the rst step the MWCS                |
 * | (Maximum Weight Common Subsequence) of the MUM sequences A and B obtained from the                | 
 * | genomes G1 and G2, where each MUM is assigned a weight (could be its length or something          |
 * | else). As in the slides, we can consider the MUM labels of one, say A, to be in canonical order   | 
 * | CS[1..n] and that of B in a permuted order PS[1..n]. An MWCS of CS[1..n] and PS[1..n] is          |
 * | in increasing order and of maximum weight.                                                        |
 * +---------------------------------------------------------------------------------------------------+
 * 
 * *****************
 * *  Test Input   *
 * *****************
 * Enter First Sequence (canonical order) A[i]: 1,2,3,4,5,6
 * Enter Second Sequence (permuted order) B[j]: 1,5,2,4,3,6
 * 
 * ****************
 * *   Output     *
 * ****************
 * 1,5,6
 */

 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
 
int N, DP[1000][1000];
char *A[1000], *B[1000];
int dire[1000][1000][3],drd=0, drl=0, drt=0;

void assignMUMWeight(void);
void computeMWCS(void);
void displayDP(void);
void showResult(void);

void parseInput(char*, char*);

/* Static Weights for test.
 * For random weights uncomment line No. 229, (that is W[i]= 1 + rnd) in  assignMUMWeight function. 
 */
int W[] = {10,2,1,3,20,2,1,3,5,4};  
 
int main(int argc, char** argv) {
    char a[1000],b[1000];
     
    printf("******************** MAXIMUM WEIGHT COMMON SUBSEQUENCE **********************\n\n");
    printf("Enter First Sequence (canonical order) A[i]: ");
    scanf("%s",&a);
    printf("\n\nEnter Second Sequence (permuted order) B[j]: ");
    scanf("%s",&b);
     
    if(strlen(a)!=strlen(b)){
        printf("\nPlease enter same length sequences as one is in canonical order and another is in permuted order.");
    }else{
        parseInput(a, b);
    }
    
    /*** Assign Random Weights to each MUMs ***/
    assignMUMWeight();
    
    
    /*** Calculate MWCS ***/
    computeMWCS();
    
    /*** Display Dynamic Programming Table ***/
    displayDP();
     
    /*** Display MWCS from DP table ***/
    showResult();
     
    return 0;
}
 
/*
 * Function: getmax()
 * Input Parameters: Three integer values 
 * Returns: Maximum among the three values.
 * And Set the flag which which is maximum
 * 
 * Here, 
 * 'd' represents the value from diagonal cell of current position(T[i][j]) in DP table, that is T[i-1][j-1].
 * 'l' represents the value from left cell of current position(T[i][j]) in DP table, that is T[i][j-1].
 * 't' represents the value from top cell of current position(T[i][j]) in DP table, that is T[i-1][j].
 * 
 * Also
 * drd, drl and drt are the global variables
 * drd=1 if the T[i-1][j-1] is maximum otherwise 0
 * drl=1 if the T[i][j-1] is maximum otherwise 0
 * drt=1 if the T[i-1][j] is maximum otherwise 0 
 */
int getmax(int d, int l, int t){
    int max=0;
    if(d==l && l==t){
        max = d;
        drd = 1;
        drl = 1;
        drt = 1;
    }else{
        if(d>l && d>t){
            max = d;
            drd = 1;
            drl = 0;
            drt = 0;
        }else if(l>d && l>t){
            max = l;
            drd = 0;
            drl = 1;
            drt = 0;
        }else{
            max = t;
            drd = 0;
            drl = 0;
            drt = 1;
        }
        
    }
    
    return max;
}

/* Function: computeMWCS
 * DP[i][j] : Store the values of Dynamic Programming table.
 * dire[i][j][] : Stores the direction from which cell(DP[i-1][j-1],DP[i-1][j],DP[i][j-1]) 
 * the current value() came.
 * 
 * dire[i][j][0]=1 if the current DP[i][j] is came from DP[i-1][j-1] otherwise 0
 * dire[i][j][1]=1 if the current DP[i][j] is came from DP[i][j-1] otherwise 0
 * dire[i][j][2]=1 if the current DP[i][j] is came from DP[i-1][j] otherwise 0
 */
void computeMWCS(){
    int i,ti,j,tj,k,d,l,t,temp,max;
     
    for(i=1; i<=N; i++){
        for(j=1; j<=N+1; j++){
            DP[i][j]=0;
        }
    }
     
    for(ti=1; ti<=N; ti++){
        for(tj=1; tj<=N; tj++){
            i=ti-1;
            j=tj-1;
            l = DP[ti][tj-1];
            t = DP[ti-1][tj];
            d = DP[ti-1][tj-1];
            if(stringToInt(A[i])==stringToInt(B[j])){  
                d += W[i];
            }
            max = getmax(d,l,t);
            DP[ti][tj] = max;
            dire[i][j][0] = drd;
            dire[i][j][1] = drl;
            dire[i][j][2] = drt;
            
        }
    }
}

/*
 * Function: showResult()
 * Description: This function find the maximum weighted path dynamic programming table 
 *              and displays the sequence 
 * 
 * dire[i][j][] : represents the direction for DP[i][j]
 * list[i] : keeps the maximum weighted common subsequence
 */

void showResult(){
    int i,j,k=0,list[1000],m=N,o=0,t;
    printf("\n\n ***********  SEQUENCE  ************\n");
    printf("Maximum WEIGHT : %d\n",DP[N][N]);
    for(i=N; i>-1; i--){
        for(j=m; j>-1; j--){
           
            if(dire[i][j][0]==1){  // if current score DP[i][j] came form diagonal cell(DP[i-1][j-1])
                list[o] = i+1;
                o++;
                m=j-1;
                break;
            }else if(dire[i][j][2]==1){  // if current score DP[i][j] came form top cell(DP[i-1][j])
                m=j;
                break;
            }else{
                
            }
        }
    }
    for(i=0; i<o;i++){
        k=o-1-i;
        if(i < k){
            t= list[i];
            list[i] = list[k];
            list[k] = t;
        }
    }
    printf("\n\n");
    for(i=0; i<o;i++){
        printf("%d ",list[i]);
    }
    
    
}

/*
 * Function: assignMUMWeight()
 * W[i] = Weight of MUM A[i] 
 * Each weight assigned randomly.
 */ 
void assignMUMWeight(){
    int i,rnd;
    for(i=0; i<N; i++){
        rnd = rand() % 20;
        //W[i]= 1 + rnd;  // assign random weight
    }
    
    printf("\n\n ----------------- MUM WEIGHTS -----------------\nMUMs         = ");
    for(i=0; i<N; i++){
        printf("%s  ",A[i]);
    }
    printf("\nWeight, W[i] = ");
    for(i=0; i<N; i++){
        printf("%d  ",W[i]);
    }
    printf("\n");
}

 /*
  * Function: displayDP()
  * This function display the Dynamic Programming table.
  * DP[i][j] = stores the score at i-th row and j-th column
  */
void displayDP(){
    int i,j;
    printf("\n\n------------------------ DP TABLE ------------------------\n");
    for(i=-1; i<=N; i++){
        if(i<0){
            char mm[]="0"; 
            printf("     ");
        }else{
            if(i==0){
                printf("0    ");  // Display first row with 0's 
            }else{
                printf("%s    ", B[i-1]);  // Display the MUMs of permuted sequence(B[i]) 
            }
             
        }
    }
    printf("\n-----------------------------------------------------------\n");
    i=0;
    for(i=-1; i<N; i++){
        if(i==-1){
            printf("0  | ");
        }else{
            printf("%s  | ",A[i]);  // Display canonical sequence(A[i])
        }
         
        for(j=0; j<=N; j++){
            printf("%d    ",DP[i+1][j]);  // Display score at DP[i][j] of Dynamic Programming Table.
        }
        printf("\n");
    }
    printf("\n\n");
}
 

/*
 * Function: parseInput()
 * This function parse the user input strings(canonical sequence and permuted sequence) 
 * and store the values in A[i](canonical seq.) and B[i](permuted seq.)
 */
void parseInput(char a[1000], char b[1000]){
    int i=0,j=0,sum=0;
     
    A[i] = strtok(a,",");
    while(A[i]!=NULL)
    {
       j= ++i;
       A[j] = strtok(NULL,",");
    }
    N=i;
    i=0;
    B[i] = strtok(b,",");
 
    while(B[i]!=NULL)
    {
       sum=0; 
       j= ++i; 
       B[j] = strtok(NULL,",");
    }
}

/*
 * Function: stringToInt()
 * Parameter: character or strings, 
 * 
 * Example: This function takes string like, '123' and 
 * converts to integer value 123 and returns it.
 */
int stringToInt(char str[]){
    int i=0,sum=0;
 
    while(str[i]!='\0'){
         if(str[i]< 48 || str[i] > 57){
             printf("Unable to convert it into integer.\n");
             return 0;
         }
         else{
             sum = sum*10 + (str[i] - 48);
             i++;
         }
    }
    return sum;
}