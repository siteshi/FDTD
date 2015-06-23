#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#include<math.h>

#define IE 60
#define JE 60

int initialize_rand_array(int ***a ,int num_row, int num_col)
{
 *a =(int **)malloc(num_row*sizeof(int *));
 int i =0,j=0;

 int *data = (int *)malloc(num_row*num_col*sizeof(int ));

  for(i=0;i<num_row*num_col;i++)
 {
	data[i] = i*i + 9 + i+2;
 }


  for(i=0;i<num_row;i++)
  {

   (*a)[i] = &data[i*num_col];
  }

 return 1;
}




int main(int argc, char *argv[])
{

 FILE *fpE,*fpH;
 MPI_Init(NULL,NULL);
 int world_rank;
 MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
 int world_size;
 MPI_Comm_size(MPI_COMM_WORLD,&world_size);
 float ga[IE][JE], dz[IE][JE], ez[IE][JE], hx[IE][JE], hy[IE][JE];
 
 int sub_IE = 60, sub_JE = 10;
 float sub_ga[ sub_IE][sub_JE], sub_dz[ sub_IE][sub_JE], sub_ez[ sub_IE][sub_JE], sub_hx[ sub_IE][sub_JE], sub_hy[ sub_IE][sub_JE];
 float recv_sub_hx[60][1], recv_sub_ez[60][1];
    int l, n, i, j, ic, jc, nsteps,  NSTEPS = 200;
    float ddx, dt, T, epsz, pi, epsilon, sigma, eaf;
    float t0, spread, pulse;
    FILE *fp = NULL;
    
    ic=IE/2;
    jc=JE/2;
    ddx=0.01;
    dt=ddx/6e8;
    epsz=8.8e-12;
    pi=3.14159;
    
   
    
    t0=30;
    spread=6.0;
    T=0;
    nsteps=1;


 /************************** this part is to intialize all global 2D array **************************************/

if(world_rank == 0)
 { 

    for (j=0; j<JE; j++){
        //printf("%2d  ",j);
        for (i=0;i<IE;i++){
            dz[i][j]=0.0;
            ez[i][j]=0.0;
            hx[i][j]=0.0;
            hy[i][j]=0.0;
            ga[i][j]=1.0;
            //printf("%5.2f ", ga[i][j]);
        }
        printf("\n");
    }



  //printf("Intial value of H and E are... \n ");
 // for(i =0;i<num;i++)
 // {
   //printf("(%2f, %2f) \n ",h_y[i],e_x[i]);
 // }
 }

/************************** End of global intialization  part **************************************/





	/************************** Define new 2D array data type **************************************/
      MPI_Datatype type, resizedtype;
      int sizes[2]    = {IE,JE};  /* size of global array */
      int subsizes[2] = {sub_IE,sub_JE};  /* size of sub-region */
      int starts[2]   = {0,0};  /* let's say we're looking at region "0",
                                 which begins at index [0,0] */

       /* as before */
       MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT , &type);
       /* change the extent of the type */
       MPI_Type_create_resized(type, 0, 10*sizeof(float ), &resizedtype);
       MPI_Type_commit(&resizedtype);

	int counts[6] = {1,1,1,1,1,1};   /* how many pieces of data everyone has, in units of blocks */
	int displs[6] = {0,1,2,3,4,5};   /* the starting point of everyone's data */
                                 /* in the global array, in block extents */
	MPI_Barrier(MPI_COMM_WORLD);

	
    	MPI_Datatype newtype;
   	int sizes_sub[2]    = {60,10};  /* size of global array */
    	int subsizes_sub[2] = {60,1};  /* size of sub-region */
    	int starts_sub[2]   = {0,0};  /* let's say we're looking at region "0",
                                 which begins at index [0,0] */

   	 MPI_Type_create_subarray(2, sizes_sub, subsizes_sub, starts_sub, MPI_ORDER_C, MPI_FLOAT, &newtype);
   	 MPI_Type_commit(&newtype);

	/************************** End of defining data type **************************************/





/************************** Main loop starts **************************************/
for (n =0 ;n<NSTEPS;n++) {

     /************************** this part is for scatter 2D array **************************************/


	MPI_Scatterv(&(ga[0][0]), counts, displs, resizedtype,&(sub_ga[0][0]), 60*10, MPI_FLOAT ,0, MPI_COMM_WORLD);

	MPI_Scatterv(&(dz[0][0]), counts, displs, resizedtype,&(sub_dz[0][0]), 60*10, MPI_FLOAT ,0, MPI_COMM_WORLD);

	MPI_Scatterv(&(ez[0][0]), counts, displs, resizedtype,&(sub_ez[0][0]), 60*10, MPI_FLOAT ,0, MPI_COMM_WORLD);

	MPI_Scatterv(&(hx[0][0]), counts, displs, resizedtype,&(sub_hx[0][0]), 60*10, MPI_FLOAT ,0, MPI_COMM_WORLD);

	MPI_Scatterv(&(hy[0][0]), counts, displs, resizedtype,&(sub_hy[0][0]), 60*10, MPI_FLOAT ,0, MPI_COMM_WORLD);


	MPI_Barrier(MPI_COMM_WORLD);
   /************************** End of scatter part **************************************/


     T=T+1;

     printf("size of world is  %d \n ",world_size);
 	for(j=0;j<sub_JE;j++){
                for(i=0;i<sub_IE;i++){
                   //printf("%5.2f ",sub_ga[i][j]); 
                }
            }
   
   /************************start of the Main loop for all the major caluclation *************************/

	if(world_rank == 0)
 	{
            printf("I am  inside processor %d and i have waiting send hx data to %d \n ",world_rank,world_rank +1);
	    MPI_Send(&(sub_hx[0][9]), 1, newtype, world_rank +1, 300, MPI_COMM_WORLD);   //MPI_Send(&(global[0][0]), 1, newtype, dest, tag, MPI_COMM_WORLD);
	    printf("I am  inside processor %d and i have succefully send hx data to %d \n ",world_rank,world_rank +1);
	    //calculate the Dz field
            for(j=1;j<sub_JE;j++){
                for(i=1;i<sub_IE;i++){
                    //// first i will receive
                    sub_dz[i][j]=sub_dz[i][j]+0.5*(sub_hy[i][j]-sub_hy[i-1][j]-sub_hx[i][j]+sub_hx[i][j-1]);
                }
            }

            //calculate the Ez field
            for(j=0;j<sub_JE;j++){
                for(i=0;i<sub_IE;i++){
                    //// first i will send  ez[1....IE][processorRank*elementPerProcessor] tag row no i.e i
                    sub_ez[i][j]=sub_ga[i][j]*sub_dz[i][j];
                }
            }
 

	   printf("I am  inside processor %d and i have waiting receive  ez data from %d \n ",world_rank,world_rank +1);
	   MPI_Recv(&(recv_sub_ez[0][0]), 60*1, MPI_FLOAT, world_rank +1, 200, MPI_COMM_WORLD,MPI_STATUS_IGNORE);  //MPI_Recv(&(local[0][0]), 3*3, MPI_INT, 0, tag, MPI_COMM_WORLD);
	   printf("I am  inside processor %d and i have succefully receive  ez data from %d \n ",world_rank,world_rank +1);
 	   for(j=0;j<1;j++){
                for(i=0;i<60;i++){
                   printf("%5.2f ",recv_sub_ez[i][j]); 
                }
            }
   
            //calculate the Hx field
            for(j=0; j<sub_JE;j++){
                for(i=0;i<sub_IE;i++){
                    //// first going to receive ez[1....IE][j+1] tag row no i.e i
		    if(j == sub_JE -1){
		    sub_hx[i][j]=sub_hx[i][j]+0.5*(sub_ez[i][j]-recv_sub_ez[i][0]);}
		    else{
		    sub_hx[i][j]=sub_hx[i][j]+0.5*(sub_ez[i][j]-sub_ez[i][j+1]);}

		   
                }
            }

            //calculate the Hy field
            for(j=0; j<sub_JE;j++){
                for(i=0;i<sub_IE-1;i++){
                    sub_hy[i][j]=sub_hy[i][j]+0.5*(sub_ez[i+1][j]-sub_ez[i][j]);
                }
            }
        

   	
	}









	if(world_rank !=0 && world_rank != world_size-1)
 	{
            printf("I am  inside processor %d and i have waiting to receive hx data from %d \n ",world_rank,world_rank -1);
	    MPI_Recv(&(recv_sub_hx[0][0]), 60*1, MPI_FLOAT, world_rank -1, 300, MPI_COMM_WORLD,MPI_STATUS_IGNORE); //MPI_Recv(&(local[0][0]), 3*3, MPI_INT, 0, tag, MPI_COMM_WORLD);
	    printf("I am  inside processor %d and i have succefully receive hx data from %d \n ",world_rank,world_rank -1);


	    printf("I am  inside processor %d and i have waiting to send hx data to %d \n ",world_rank,world_rank +1);
            MPI_Send(&(sub_hx[0][9]), 1, newtype, world_rank +1, 300, MPI_COMM_WORLD);   //MPI_Send(&(global[0][0]), 1, newtype, dest, tag, MPI_COMM_WORLD);
	    printf("I am  inside processor %d and i have succefully send hx data to %d \n ",world_rank,world_rank +1);


	    //calculate the Dz field
            for(j=0;j<sub_JE;j++){
                for(i=1;i<sub_IE;i++){
                    //// first i will receive
		    if(j ==0 ){
                    sub_dz[i][j]=sub_dz[i][j]+0.5*(sub_hy[i][j]-sub_hy[i-1][j]-sub_hx[i][j]+recv_sub_hx[i][0]); }
		    else{
                    sub_dz[i][j]=sub_dz[i][j]+0.5*(sub_hy[i][j]-sub_hy[i-1][j]-sub_hx[i][j]+sub_hx[i][j-1]); }

                    
                }
            }



           //put a Guassian pulse in the middle of global array
	   if(world_rank == world_size/2 ){
            
            pulse=exp(-0.5*pow( (t0-T)/spread,2.0 ) );
            sub_dz[sub_IE/2][2]=pulse;
	    }



            //calculate the Ez field
            for(j=0;j<sub_JE;j++){
                for(i=0;i<sub_IE;i++){
                    //// first i will send  ez[1....IE][processorRank*elementPerProcessor] tag row no i.e i
                    sub_ez[i][j]=sub_ga[i][j]*sub_dz[i][j];
                }
            }


           printf("I am  inside processor %d and i have waiting receive ez data from %d \n ",world_rank,world_rank +1);
	   MPI_Recv(&(recv_sub_ez[0][0]), 60*1, MPI_FLOAT, world_rank +1, 200, MPI_COMM_WORLD,MPI_STATUS_IGNORE);  //MPI_Recv(&(local[0][0]), 3*3, MPI_INT, 0, tag, MPI_COMM_WORLD);
	   printf("I am  inside processor %d and i have succefully receive ez data from %d \n ",world_rank,world_rank +1);
	   printf("I am  inside processor %d and i have waiting send ez data to %d \n ",world_rank,world_rank -1);
	   MPI_Send(&(sub_ez[0][0]), 1, newtype, world_rank -1, 200, MPI_COMM_WORLD);   //MPI_Send(&(global[0][0]), 1, newtype, dest, tag, MPI_COMM_WORLD);
	   printf("I am  inside processor %d and i have succefully send ez data to %d \n ",world_rank,world_rank -1);

	  	   

            //calculate the Hx field
            for(j=0; j<sub_JE;j++){
                for(i=0;i<sub_IE;i++){
                    //// first going to receive ez[1....IE][j+1] tag row no i.e i
		    if(j == sub_JE -1){
                     sub_hx[i][j]=sub_hx[i][j]+0.5*(sub_ez[i][j]-recv_sub_ez[i][0]);}
		    else{
                     sub_hx[i][j]=sub_hx[i][j]+0.5*(sub_ez[i][j]-sub_ez[i][j+1]);}

		   
                }
            }

            //calculate the Hy field
            for(j=0; j<sub_JE;j++){
                for(i=0;i<sub_IE-1;i++){
                    sub_hy[i][j]=sub_hy[i][j]+0.5*(sub_ez[i+1][j]-sub_ez[i][j]);
                }
            }

   	
	}



	if(world_rank == world_size-1)
 	{
            printf("I am  inside processor %d and i have waiting to receive hx data from %d \n ",world_rank,world_rank -1);
	    MPI_Recv(&(recv_sub_hx[0][0]), 60*1, MPI_FLOAT, world_rank -1, 300, MPI_COMM_WORLD,MPI_STATUS_IGNORE); //MPI_Recv(&(local[0][0]), 3*3, MPI_INT, 0, tag, MPI_COMM_WORLD);
            printf("I am  inside processor %d and i have receive hx data from %d \n ",world_rank,world_rank -1);


	    //calculate the Dz field
            for(j=0;j<sub_JE;j++){
                for(i=1;i<sub_IE;i++){
                    //// first i will receive
		    if(j ==0 ){
                    sub_dz[i][j]=sub_dz[i][j]+0.5*(sub_hy[i][j]-sub_hy[i-1][j]-sub_hx[i][j]+recv_sub_hx[i][0]); }
		    else{
                    sub_dz[i][j]=sub_dz[i][j]+0.5*(sub_hy[i][j]-sub_hy[i-1][j]-sub_hx[i][j]+sub_hx[i][j-1]); }

                    
                }
            }

            //calculate the Ez field
            for(j=0;j<sub_JE;j++){
                for(i=0;i<sub_IE;i++){
                    //// first i will send  ez[1....IE][processorRank*elementPerProcessor] tag row no i.e i
                    sub_ez[i][j]=sub_ga[i][j]*sub_dz[i][j];
			
                }
            }

	   printf("I am  inside processor %d and i have waiting send ez data to %d \n ",world_rank,world_rank -1);
	   MPI_Send(&(sub_ez[0][0]), 1, newtype, world_rank -1, 200, MPI_COMM_WORLD); //MPI_Send(&(global[0][0]), 1, newtype, dest, tag, MPI_COMM_WORLD);
	   printf("I am  inside processor %d and i have succefully send ez data to %d \n ",world_rank,world_rank -1);
            //calculate the Hx field
            for(j=0; j<sub_JE;j++){
                for(i=0;i<sub_IE;i++){
                    //// first going to receive ez[1....IE][j+1] tag row no i.e i
		    if(j == sub_JE -1){
                     sub_hx[i][j]=sub_hx[i][j]+0.5*(sub_ez[i][j]-0);}   //sub_hx[i][j]=sub_hx[i][j]+0.5*(sub_ez[i][j]-recv_sub_ez[i][0]);
		    else{
                     sub_hx[i][j]=sub_hx[i][j]+0.5*(sub_ez[i][j]-sub_ez[i][j+1]);}

		   
                }
            }

            //calculate the Hy field
            for(j=0; j<sub_JE;j++){
                for(i=0;i<sub_IE-1;i++){
                    sub_hy[i][j]=sub_hy[i][j]+0.5*(sub_ez[i+1][j]-sub_ez[i][j]);
                }
            }

   	
	}



MPI_Barrier(MPI_COMM_WORLD);

     /************************** this part is for gather 2D array **************************************/

        MPI_Gatherv(&(sub_ga[0][0]), 60*10, MPI_FLOAT,&(ga[0][0]),counts,displs, resizedtype,0, MPI_COMM_WORLD);
	MPI_Gatherv(&(sub_dz[0][0]), 60*10, MPI_FLOAT,&(dz[0][0]),counts,displs, resizedtype,0, MPI_COMM_WORLD);
        MPI_Gatherv(&(sub_ez[0][0]), 60*10, MPI_FLOAT,&(ez[0][0]),counts,displs, resizedtype,0, MPI_COMM_WORLD);
	MPI_Gatherv(&(sub_hx[0][0]), 60*10, MPI_FLOAT,&(hx[0][0]),counts,displs, resizedtype,0, MPI_COMM_WORLD);
	MPI_Gatherv(&(sub_hy[0][0]), 60*10, MPI_FLOAT,&(hy[0][0]),counts,displs, resizedtype,0, MPI_COMM_WORLD);

   /************************** End of gather part **************************************/
	if(world_rank ==0)
	{
            
  
  	char fName[20] = "Ez";
  	char str[10];
   	sprintf(str, "%d", n);
   	strcat(fName,str);
            
   	fp=fopen(fName,"w");
            for(i=0;i<IE;i++){
                for(j=0;j<JE;j++){
                    fprintf(fp, "%6.3f,",ez[i][j]);
                }
                fprintf(fp,"\n");
            }
            
            //fprintf(fp,"\n");
            fclose(fp);
	}


}
/************************** main loop ends **************************************/
 
  MPI_Finalize();

   }

















