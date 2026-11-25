#include "iostream"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "fstream"
#include "vector"
#include "string"
using namespace std;




#define M 746*3                                // number of variables
#define N 746*5+1 
//int N;
//int M;
#define matrice "autotroph_matrix_configuration_1.dat"
#define minimum 1e-9
std::vector<int> matrix[M];
std::vector<long double > matrix2[M];
 
void reading(){
   fstream file;
   file.open(matrice, ios::in);                                    
   int  i1=0;
   int  j1=0;        
   long double  pepe;
   while(file >> pepe){
          if(pepe!=0){
                     matrix[i1].push_back(j1);
                     matrix2[i1].push_back(pepe);
                     }
          j1++;
          if(j1==N){
                   j1=0;
                   i1++;
                   }
          }        
   file.close();     
   }

void subtractpiv( int k1, int k2){
                                long double  riga[N];
                                for(int i=0;i<=N-1;i++) riga[i]=0;
                                long double   g = matrix2[k2][0]/matrix2[k1][0];
                                for(int i=1;i<=matrix[k1].size()-1;i++) riga[matrix[k1][i]]  = matrix2[k1][i]*g;
                                for(int i=1;i<=matrix[k2].size()-1;i++) riga[matrix[k2][i]] -= matrix2[k2][i];
                                matrix[k1].clear();
                                matrix2[k1].clear();
                                for(int i=0;i<=N-1;i++){
                                                  if(fabs(riga[i])>minimum){
                                                        matrix[k1].push_back(i);
                                                        matrix2[k1].push_back(riga[i]);
                                                        }
                                                  }
                              }

   



  
void quicksort(int k, int km, int orders[], int pivots[]){
                            int centre;
                            if(k-km>=1){
                               int pivot = km+int((k-km)/2);
                               int l,p;
                               l=0;
                               p=k-km-1;
                               long double  neworders[k-km];
                               for(int i=km;i<=k-1;i++){
                                                   if(i != pivot){
                                                                      if(pivots[orders[i]]<pivots[orders[pivot]]){
                                                                               neworders[l]=orders[i];
                                                                               l++;
                                                                              }
                                                                      else{
                                                                               neworders[p]=orders[i];
                                                                               p--;
                                                                           }
                                                              }
                                                   }
                               
                               neworders[p]=orders[pivot];
                               for(int i=km;i<=k-1;i++) orders[i]=neworders[i-km];
                               centre = p+km;
                               quicksort(k,centre+1,orders,pivots);
                               quicksort(centre,km, orders,pivots);
                               }   
                                             
                               
                               
   
                               }
         





int main (int argc, char* argv[]){
     cout.precision(6);

     //M=atof(int argv[1]);
     //N=atof(int argv[2]);
     //cout << M;
     //int M=argv[1];
     //int N=argv[2];
     
     //std::vector<int> matrix[M];
     //std::vector<long double> matrix2[M];
 
 
     reading();
     

     int ok=0;
     int order[M];
     for(int i=0;i<=M-1;i++) order[i]=i;
     int pivots[M];
     for(int i=0;i<=M-1;i++){
                 if(matrix[i].size()>0) pivots[i]=matrix[i][0];
                 else pivots[i]=80000;
                                
                 }
     
     while(ok==0){
        
       quicksort(M,0,order,pivots);
       for(int j=0;j<M-1;j++){
         if(pivots[order[j+1]]==pivots[order[j]] && pivots[order[j]]!=80000){ 
          long double min1=100000000;
          if(matrix[order[j]].size()>1) for(int i=0;i<=matrix[order[j]].size()-1;i++) if(fabs(matrix2[order[j]][0]/matrix2[order[j]][i])<min1) min1 = fabs(matrix2[order[j]][0]/matrix2[order[j]][i]);
          long double min2=100000000;
          if(matrix[order[j+1]].size()>1) for(int i=0;i<=matrix[order[j+1]].size()-1;i++)if(fabs(matrix2[order[j+1]][0]/matrix2[order[j+1]][i])<min2) min2 = fabs(matrix2[order[j+1]][0]/matrix2[order[j+1]][i]); 
          if(min2>min1){
                 int k2 = order[j+1];
                 order[j+1]=order[j];
                 order[j] = k2;              
                 }
           }
         }
       ok=1;
       int j=0;
       for(int j=0;j<M-1;j++){

                      if(pivots[order[j+1]]==pivots[order[j]] && pivots[order[j]]!=80000){
                                        
                                        subtractpiv(order[j+1],order[j]);           
                                        ok=0;
                                        if(matrix[order[j+1]].size()>0) pivots[order[j+1]]=matrix[order[j+1]][0];
                                        else pivots[order[j+1]]=80000;
                                        }
                          
                             }
        
       };
     
    for(int i=0;i<=M-1;i++) {
                            if(matrix[i].size()>0){
                                      long double  norm = matrix2[i][0];
                                      for(int j=0;j<=matrix[i].size()-1;j++) matrix2[i][j]/=norm;
                                      }
                          } 


    

    
    for(int k1=M-2;k1>=0;k1--){
                           
                           int k = order[k1];
                           if(matrix[k].size()>1){
                              for(int i=1;i<=matrix[k].size()-1;i++){
                                               for(int j1=k1+1;j1<=M-1;j1++){
                                                                int j = order[j1];
                                                                if(matrix[j].size()>0){
                                                                    if(matrix[j][0]==matrix[k][i]){
                                                                                              long double  rigak[N];
                                                                                              for(int a=0;a<=N-1;a++) rigak[a]=0;
                                                                                              for(int a=0;a<=matrix[k].size()-1;a++) rigak[matrix[k][a]] = matrix2[k][a];        
                                                                                              for(int a=0;a<=matrix[j].size()-1;a++) rigak[matrix[j][a]]-= matrix2[j][a]*matrix2[k][i];
                                                                                              matrix[k].clear();
                                                                                              matrix2[k].clear(); 
                                                                                              for(int a=0;a<=N-1;a++){
                                                                                                            if(rigak[a]!=0){
                                                                                                                   matrix[k].push_back(a); 
                                                                                                                   matrix2[k].push_back(rigak[a]);
                                                                                                                   }
                                                                                                            }
                                                                                               }
                                                                                       }
                                                                       }           
                                                                   }
                                                 }
                                 }
       
     int indip[N];
     for(int i=0;i<=N-1;i++) indip[i]=M+1;
     for(int i=0;i<=M-1;i++) if(matrix[i].size()>0) indip[matrix[i][0]]=i;
     int q=0;
     for(int i=0;i<=N-1;i++) if(indip[i]==M+1) q++;
     int quelo=0;
     for(int k=0;k<=N-1;k++){
         if(indip[k]!=M+1){
              long double  riga[N];   
              int k1 = indip[k];
              for(int j=0;j<=N-1;j++) riga[j]=0;
              for(int j=0;j<=matrix[k1].size()-1;j++) riga[matrix[k1][j]] = matrix2[k1][j]; 
              for(int i=0;i<=N-1;i++) if(indip[i]==M+1){
                                                            if(fabs(riga[i])>minimum){
                                                                                    int s=1;
                                                                                    if(riga[i]<0) s=-1;
                                                                                    if(fabs(riga[i]-int(riga[i]+s*minimum))<minimum) riga[i] = int(riga[i]+s*minimum);
                                                                                    
                                                                                    cout << -riga[i] << "  ";
                                                                                    }
                                                            else cout << 0 << "  ";
                                                       }
              cout << endl;                                    
              }
         else{
             for(int j=0;j<=q-1;j++){
                               int c=0;
                               if(quelo==j) c=1; 
                               cout << c << "  ";
                               }
              quelo++;
              cout << endl;
            }
       
        }
      



      return 0;
 }
