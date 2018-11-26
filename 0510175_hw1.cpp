#include <iostream>
#include <ctime>
#include "fstream"
using namespace std;



template<typename T>
class Strassen_class{
public:
    
    void ADD(T** MatrixA, T** MatrixB, T** MatrixResult, int MatrixSize );
    void SUB(T** MatrixA, T** MatrixB, T** MatrixResult, int MatrixSize );
    void MUL( T** MatrixA, T** MatrixB, T** MatrixResult, int MatrixSize ); //暴力演算法
    void FillMatrix( T** MatrixA, T** MatrixB,int argc, char *argv[]);                 //A,B矩陣給值
    void PrintMatrix(T **MatrixA,int MatrixSize);                           //印出矩陣
    void Strassen(int N, T **MatrixA, T **MatrixB, T **MatrixC);            //Strassen演算法
   
    
    void Expend(T** MatrixA, T** MatrixB, int oldsize,int newsize);                //矩陣擴增if n!=2乘冪
    //void Expend2(T** MatrixA, T** MatrixB, int MatrixSize);
    
    
};





inline bool is_pow_of_two(int i)
{
    if (i < 1) return false;
    return (i & (i - 1)) == 0;
}





int reshape(int MatrixSize){
    while(is_pow_of_two(MatrixSize)!=1){
        MatrixSize=++MatrixSize;
    }


    return MatrixSize ;
}




template<typename T>
void Strassen_class<T>::Expend(T** MatrixA, T** MatrixB, int oldsize,int newsize){

    
    
    for (int i=oldsize; i<=newsize; i++) {
        for (int j=0; j<oldsize; j++) {
           MatrixB[i][j] = (MatrixA[i][j] = 0);
        }
    }
    
    for (int j=oldsize; j<=newsize; j++) {
        for (int i=0; i<=oldsize; i++) {
            MatrixB[i][j] = (MatrixA[i][j] = 0);
            
        }
    }
    
    ////////擺右下角單位矩陣
    for (int i=oldsize; i<=newsize; i++) {
        
        for (int j=oldsize; j<=newsize; j++) {
            if(i==j){
                MatrixB[i][j] = (MatrixA[i][j] = 1);
                //cout<<i<<endl<<j<<endl;
            }
            else{
                MatrixB[i][j] = (MatrixA[i][j] = 0);
                //cout<<i<<endl;
            }
        }
    }
    /////擺右下角單位矩陣
    
    }



template<typename T>
void Strassen_class<T>::ADD(T** MatrixA, T** MatrixB, T** MatrixResult, int MatrixSize )
{
    for ( int i = 0; i < MatrixSize; i++)
    {
        for ( int j = 0; j < MatrixSize; j++)
        {
            MatrixResult[i][j] =  MatrixA[i][j] + MatrixB[i][j];
        }
    }
}





template<typename T>
void Strassen_class<T>::SUB(T** MatrixA, T** MatrixB, T** MatrixResult, int MatrixSize )
{
    for ( int i = 0; i < MatrixSize; i++)
    {
        for ( int j = 0; j < MatrixSize; j++)
        {
            MatrixResult[i][j] =  MatrixA[i][j] - MatrixB[i][j];
        }
    }
}





template<typename T>
void Strassen_class<T>::MUL( T** MatrixA, T** MatrixB, T** MatrixResult, int MatrixSize )
{
    for (int i=0;i<MatrixSize ;i++)
    {
        for (int j=0;j<MatrixSize ;j++)
        {
            MatrixResult[i][j]=0;
            for (int k=0;k<MatrixSize ;k++)
            {
                MatrixResult[i][j]=MatrixResult[i][j]+MatrixA[i][k]*MatrixB[k][j];
            }
        }
    }
}






template<typename T>
void Strassen_class<T>::Strassen(int N, T **MatrixA, T **MatrixB, T **MatrixC)
{
    
    int HalfSize = N/2;
    int newSize = N/2;
    
    if ( N <= 64 )    //小於64用暴力演算法
    {
        MUL(MatrixA,MatrixB,MatrixC,N);
    }
    else
    {
        T** A11;
        T** A12;
        T** A21;
        T** A22;
        
        T** B11;
        T** B12;
        T** B21;
        T** B22;
        
        T** C11;
        T** C12;
        T** C21;
        T** C22;
        
        T** M1;
        T** M2;
        T** M3;
        T** M4;
        T** M5;
        T** M6;
        T** M7;
        T** AResult;
        T** BResult;
        
        //making a 1 diminsional pointer based array.
        A11 = new T *[newSize];
        A12 = new T *[newSize];
        A21 = new T *[newSize];
        A22 = new T *[newSize];
        
        B11 = new T *[newSize];
        B12 = new T *[newSize];
        B21 = new T *[newSize];
        B22 = new T *[newSize];
        
        C11 = new T *[newSize];
        C12 = new T *[newSize];
        C21 = new T *[newSize];
        C22 = new T *[newSize];
        
        M1 = new T *[newSize];
        M2 = new T *[newSize];
        M3 = new T *[newSize];
        M4 = new T *[newSize];
        M5 = new T *[newSize];
        M6 = new T *[newSize];
        M7 = new T *[newSize];
        
        AResult = new T *[newSize];
        BResult = new T *[newSize];
        
        int newLength = newSize;
        
        //making that 1 diminsional pointer based array , a 2D pointer based array
        for ( int i = 0; i < newSize; i++)
        {
            A11[i] = new T[newLength];
            A12[i] = new T[newLength];
            A21[i] = new T[newLength];
            A22[i] = new T[newLength];
            
            B11[i] = new T[newLength];
            B12[i] = new T[newLength];
            B21[i] = new T[newLength];
            B22[i] = new T[newLength];
            
            C11[i] = new T[newLength];
            C12[i] = new T[newLength];
            C21[i] = new T[newLength];
            C22[i] = new T[newLength];
            
            M1[i] = new T[newLength];
            M2[i] = new T[newLength];
            M3[i] = new T[newLength];
            M4[i] = new T[newLength];
            M5[i] = new T[newLength];
            M6[i] = new T[newLength];
            M7[i] = new T[newLength];
            
            AResult[i] = new T[newLength];
            BResult[i] = new T[newLength];
            
            
        }
        //splitting input Matrixes, into 4 submatrices each.
        for (int i = 0; i < N / 2; i++)
        {
            for (int j = 0; j < N / 2; j++)
            {
                A11[i][j] = MatrixA[i][j];
                A12[i][j] = MatrixA[i][j + N / 2];
                A21[i][j] = MatrixA[i + N / 2][j];
                A22[i][j] = MatrixA[i + N / 2][j + N / 2];
                
                B11[i][j] = MatrixB[i][j];
                B12[i][j] = MatrixB[i][j + N / 2];
                B21[i][j] = MatrixB[i + N / 2][j];
                B22[i][j] = MatrixB[i + N / 2][j + N / 2];
                
            }
        }
        
        //here we calculate M1..M7 matrices .
        //M1[][]
        ADD( A11,A22,AResult, HalfSize);
        ADD( B11,B22,BResult, HalfSize);                //p5=(a+d)*(e+h)
        Strassen( HalfSize, AResult, BResult, M1 ); //now that we need to multiply this , we use the strassen itself .
        
        
        //M2[][]
        ADD( A21,A22,AResult, HalfSize);              //M2=(A21+A22)B11   p3=(c+d)*e
        Strassen(HalfSize, AResult, B11, M2);       //Mul(AResult,B11,M2);
        
        //M3[][]
        SUB( B12,B22,BResult, HalfSize);              //M3=A11(B12-B22)   p1=a*(f-h)
        Strassen(HalfSize, A11, BResult, M3);       //Mul(A11,BResult,M3);
        
        //M4[][]
        SUB( B21, B11, BResult, HalfSize);           //M4=A22(B21-B11)    p4=d*(g-e)
        Strassen(HalfSize, A22, BResult, M4);       //Mul(A22,BResult,M4);
        
        //M5[][]
        ADD( A11, A12, AResult, HalfSize);           //M5=(A11+A12)B22   p2=(a+b)*h
        Strassen(HalfSize, AResult, B22, M5);       //Mul(AResult,B22,M5);
        
        
        //M6[][]
        SUB( A21, A11, AResult, HalfSize);
        ADD( B11, B12, BResult, HalfSize);             //M6=(A21-A11)(B11+B12)   p7=(c-a)(e+f)
        Strassen( HalfSize, AResult, BResult, M6);    //Mul(AResult,BResult,M6);
        
        //M7[][]
        SUB(A12, A22, AResult, HalfSize);
        ADD(B21, B22, BResult, HalfSize);             //M7=(A12-A22)(B21+B22)    p6=(b-d)*(g+h)
        Strassen(HalfSize, AResult, BResult, M7);     //Mul(AResult,BResult,M7);
        
        //C11 = M1 + M4 - M5 + M7;
        ADD( M1, M4, AResult, HalfSize);
        SUB( M7, M5, BResult, HalfSize);
        ADD( AResult, BResult, C11, HalfSize);
        
        //C12 = M3 + M5;
        ADD( M3, M5, C12, HalfSize);
        
        //C21 = M2 + M4;
        ADD( M2, M4, C21, HalfSize);
        
        //C22 = M1 + M3 - M2 + M6;
        ADD( M1, M3, AResult, HalfSize);
        SUB( M6, M2, BResult, HalfSize);
        ADD( AResult, BResult, C22, HalfSize);
        
        
        //組合小矩陣到大矩陣
        for (int i = 0; i < N/2 ; i++)
        {
            for (int j = 0 ; j < N/2 ; j++)
            {
                MatrixC[i][j] = C11[i][j];
                MatrixC[i][j + N / 2] = C12[i][j];
                MatrixC[i + N / 2][j] = C21[i][j];
                MatrixC[i + N / 2][j + N / 2] = C22[i][j];
            }
        }
        
        // 釋放矩陣記憶體空間
        for (int i = 0; i < newLength; i++)
        {
            delete[] A11[i];delete[] A12[i];delete[] A21[i];
            delete[] A22[i];
            
            delete[] B11[i];delete[] B12[i];delete[] B21[i];
            delete[] B22[i];
            delete[] C11[i];delete[] C12[i];delete[] C21[i];
            delete[] C22[i];
            delete[] M1[i];delete[] M2[i];delete[] M3[i];delete[] M4[i];
            delete[] M5[i];delete[] M6[i];delete[] M7[i];
            delete[] AResult[i];delete[] BResult[i] ;
        }
        delete[] A11;delete[] A12;delete[] A21;delete[] A22;
        delete[] B11;delete[] B12;delete[] B21;delete[] B22;
        delete[] C11;delete[] C12;delete[] C21;delete[] C22;
        delete[] M1;delete[] M2;delete[] M3;delete[] M4;delete[] M5;
        delete[] M6;delete[] M7;
        delete[] AResult;
        delete[] BResult ;
        
    }//end of else
    
}







template<typename T>
void Strassen_class<T>::PrintMatrix(T **MatrixA,int MatrixSize)
{
    ofstream output("output.txt",ios::out);
    output<<endl;
    
    for(int row = 0; row<MatrixSize; row++)
    {
        for(int column = 0; column<MatrixSize; column++)
        {
            
            
            output<<MatrixA[row][column]<<"\t";
            if ((column+1)%((MatrixSize)) == 0)
                output<<endl;
        }
        
    }
//    output.close();
    output<<endl;
}

template<typename T>
void Strassen_class<T>::FillMatrix( T** MatrixA, T** MatrixB,int argc, char *argv[])
{
    ifstream input(argv[1],ios::in);      ////測試用case0寫死 要改回input(argv[1],ios::in),fill and main
    //ifstream input("case0.txt",ios::in);
    int N;
    input >> N;
    
    ////////input
    for(int row = 0; row<N; row++)
    {
        for(int column = 0; column<N; column++)
        {
            input >> MatrixA[row][column];
        }
    }
    
    for(int row = 0; row<N; row++)
    {
        for(int column = 0; column<N; column++)
        {
            input >> MatrixB[row][column];
        }
    }
    
    input.close();
    /////////
    
}



int main(int argc, char *argv[])
{
    Strassen_class<int> stra;           //宣告Strassen_class類別
    
    
    ifstream input(argv[1],ios::in);    //測試用血死case0 要改回input(argv[1],ios::in),fill and main
    //ifstream input("case0.txt",ios::in);
    ofstream output("output.txt",ios::out);
    int N;
    input >> N;
    
    int MatrixSize = N;                     ///為可讀性
    int newSize=reshape(MatrixSize);        ///改寫矩陣維度直到2乘冪(太浪費記憶體)
    
    
    int** MatrixA;                      //宣告矩陣A
    int** MatrixB;                      //宣告矩陣B
    int** MatrixC;                      //宣告结果
    
    
    
    //////初始化pointer array
    MatrixA = new int *[MatrixSize];
    MatrixB = new int *[MatrixSize];
    MatrixC = new int *[MatrixSize];
    
    
    
    for (int i = 0; i < newSize; i++)
    {
        MatrixA[i] = new int [MatrixSize];
        MatrixB[i] = new int [MatrixSize];
        MatrixC[i] = new int [MatrixSize];
    }
    /////////
    
    
    stra.FillMatrix(MatrixA,MatrixB,argc,argv);               //矩陣給值(pass by reference)
    //stra.Expend(MatrixA, MatrixB, MatrixSize,newSize);      //擴充矩陣到2的乘冪次
    
    
    
//    clock_t vio_startTime=0 ;
//    clock_t vio_endTime =0;
//    clock_t vio_duration= vio_endTime-vio_startTime;
//
//    clock_t stra_startTime = 0 ;
//    clock_t stra_endTime = 0 ;
//    clock_t stra_duration=stra_endTime-stra_startTime;
//    srand(time_t(0));
    
    

    //*******************Strassen 乘法
    if(is_pow_of_two(N))
    {
        //cout<<"\nStrassen演算法開始時間: "<< (stra_startTime = clock());
        stra.Strassen(MatrixSize, MatrixA, MatrixB, MatrixC ); //strassen矩陣相乘算法
        //cout<<"\nStrassen演算法結束時間: "<<(stra_endTime = clock());
        //cout<<"\n運算結果... \n";
        stra.PrintMatrix(MatrixC,MatrixSize);
        
    }else{
        
        //*******************暴力演算法
        //cout<<"\n暴力演算法開始時間:  "<< (vio_startTime = clock());
        stra.MUL(MatrixA,MatrixB,MatrixC,MatrixSize);//暴力矩陣相乘算法 T(n) = O(n^3)
        //cout<<"\n暴力演算法結束時間: "<< (vio_endTime = clock());
        //cout<<"\n運算結果... \n";
        
        
        stra.PrintMatrix(MatrixC,MatrixSize);
    }

    


    //cout<<"矩陣大小 "<<MatrixSize;
    //cout<<"\n暴力矩陣演算法: "<<(vio_duration)<<" Clocks.."<<(vio_duration)/CLOCKS_PER_SEC<<" Sec";
    //cout<<"\nStrassen演算法:"<<(stra_duration)<<" Clocks.."<<(stra_duration)/CLOCKS_PER_SEC<<" Sec\n";

    
    
    
    return 0;
}

