
#include "luheader.cuh"

void linearfind(double *ans, double *matL, double *matU,double *matb,double *tempd, int size){
  int i,j;
  double d,x;
  
  //forward substitution
  for (i=0;i<size;i++){
    d=matb[i];
    for (j=0;j<i;j++){
      d=d-matL[i*size+j]*tempd[j]; 
    }
    tempd[i]=d;
  }
  
  //backward substitution
  for (i=size-1;i>=0;i--){
    x=tempd[i]; 
    for (j=i+1;j<size;j++){
      x=x-matU[i*size+j]*ans[j]; 
    }
    ans[i]=x/matU[i*size+i];
  }	
  
}

int main(int argc, char *argv[]){
  int size;
  
  if(argc<4){
    fprintf(stderr,"Please enter all args eg: matA.txt ans.txt report.txt \n where the format is Ax=b ");
    exit(EXIT_FAILURE);
    }
    
  //reading files

  std::ifstream file(argv[1]);
  if (!file.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    exit(EXIT_FAILURE);}

  file >> size;
  
  double *matA=(double*)malloc(size*size*sizeof(double));
  double *matL=(double*)malloc(size*size*sizeof(double));
  double *matU=(double*)malloc(size*size*sizeof(double));
  double *matb=(double*)malloc(size*sizeof(double));
  double *tempd=(double*)malloc(size*sizeof(double));
  double *ans=(double*)malloc(size*sizeof(double));

  //read start for matrix A & B
  clock_t ra=clock(); 

  for (int i = 0; i < size * size; i++) {
      file >> matA[i];
  }
  for (int i = 0; i < size; i++) {
      file >> matb[i];
  }

  file.close();

  //time for read A &B
  clock_t re=clock(); //read end for matrix
  double rt=(double)((re-ra)/(double)CLOCKS_PER_SEC);

  // declaration of output file
  FILE *fp1;
  fp1=fopen(argv[3],"w");
  fprintf(fp1,"Time taken to read A and B matrices from file: %1.10f s\n",rt);

  

  int i;

  
  // linear equation total time with start and stop
  clock_t start=clock();

  
  ludecompose(matL,matU,matA,size,fp1);// in gpu
  
  linearfind(ans, matL, matU,matb,tempd,size); //in cpu
  
  clock_t stop=clock();
    double cputime=(double)((stop-start)/(double)CLOCKS_PER_SEC);
    fprintf(fp1,"Total time taken in solving system of linear equations for the given size: %1.10f s\n",cputime);
  
    // printarray(ans,size,1);

  //writing to ans.txt file
  FILE *fp2;
  fp2=fopen(argv[2],"w");
  // fprintf(fp2,"size:\n");
  fprintf(fp2,"%d\n",size);

  // fprintf(fp2,"\nLower triangular matrix:\n\n");
  for(i=0;i<(size);i++)
  {
    for(int j=0;j<size;j++)
    {
      fprintf(fp2,"%f ",matL[i*size+j]);
    }
    fprintf(fp2,"\n");
  }

  // fprintf(fp2,"\nUpper triangular matrix:\n\n");
  fprintf(fp2,"\n");
  for(i=0;i<(size);i++)
  {
    for(int j=0;j<size;j++)
    {
      fprintf(fp2,"%f ",matU[i*size+j]);
    }
    fprintf(fp2,"\n");
  }

  // fprintf(fp2,"\nSolution matrix:\n\n");
  for (i=0;i<size;i++){
    fprintf(fp2,"%.10f ",ans[i]);
    fprintf(fp2,"\n");
  }	
  fclose(fp2);	
  fclose(fp1);
  
  
  return 0;

}