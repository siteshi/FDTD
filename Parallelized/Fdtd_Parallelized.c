#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#include<math.h>

float *initialize_array(int num_element)
{
 float* num_array =(float*)malloc(num_element*sizeof(float));
 int i =0;

 for(i =0;i< num_element;i++)
  {
   //num_array[i] = rand()/(float)RAND_MAX;
     num_array[i] = 0;
  }

 return num_array;
}





int main(int argc, char *argv[])
{

 
 MPI_Init(NULL,NULL);
 int world_rank;
 MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
 int world_size;
 MPI_Comm_size(MPI_COMM_WORLD,&world_size);
 int num = 200;
 int num_element_per_proc = num/world_size;
 int i,j;

 FILE *fpE,*fpH;
 

 float *h_y= NULL;
 float *e_x= NULL;
 float *final_e = NULL;
 float *final_h = NULL;

 int n, k, kc, ke, NSTEPS = 500,KE =num_element_per_proc ;
 
 float t0, spread, pulse;
 float T;
 t0 = 40;
 spread = 12;
 float data_send,data_recv;
 kc = KE/2;

 if(world_rank == 0)
 { 
  //fpE = fopen( "E-field", "w");
  fpH = fopen( "H-field", "w");
  //fpH = fopen( "Hy", "w");
  h_y = initialize_array(num);
  e_x = initialize_array(num);
  printf("Intial value of H and E are... \n ");
  for(i =0;i<num;i++)
  {
   printf("(%2f, %2f) \n ",h_y[i],e_x[i]);
  }
 }

 float *sub_h_y= (float*)malloc(num_element_per_proc*sizeof(float));
 float *sub_e_x= (float*)malloc(num_element_per_proc*sizeof(float));

 for (n =0 ;n<NSTEPS;n++) {

 MPI_Scatter(h_y,num_element_per_proc,MPI_FLOAT,sub_h_y,num_element_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
 MPI_Scatter(e_x,num_element_per_proc,MPI_FLOAT,sub_e_x,num_element_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);

 MPI_Barrier(MPI_COMM_WORLD); // this will insure that both scatter are finsihed over each slaves

 //printf("data inside processor %d \n ",world_rank);

 for(i=0;i<num_element_per_proc;i++)
 {
  printf("data inside processor %d is (%2f, %2f) \n ",world_rank,sub_h_y[i],sub_e_x[i]);
 }


 MPI_Barrier(MPI_COMM_WORLD); // this is just for test 



  if(world_rank == 0)
 {
   
   sub_e_x[0] = sub_e_x[0] + 0.5*(0.0 - sub_h_y[0]);
   data_send = sub_h_y[KE-1];
   MPI_Send(&data_send,1,MPI_FLOAT,world_rank +1 ,100,MPI_COMM_WORLD);    // 100 tag for h_x
   printf("I am  inside processor %d and i have succefully send data to %d and data is %f \n ",world_rank,world_rank +1,data_send);
   	for ( k = 1; k < KE; k++ ) {
          sub_e_x[k] = sub_e_x[k] + 0.5*(sub_h_y[k-1] - sub_h_y[k]);
        }
        //printf("I am  inside processor %d  and i have succefully calculate sub_ex \n",world_rank);
   
	MPI_Recv(&data_recv,1,MPI_FLOAT, world_rank + 1,200,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // 200 tag for e_x
        printf("I am  inside processor %d and i have succefully recv data from %d and data is %f\n ",world_rank,world_rank +1,data_recv);

        /* Calculate the Hy field */
        for ( k = 0; k < KE-1; k++ ) {
          sub_h_y[k] = sub_h_y[k] + .5*(sub_e_x[k] - sub_e_x[k+1] );
        }
	sub_h_y[k] = sub_h_y[k] + 0.5*( sub_e_x[k]- data_recv );
        
   // printf("I am  inside processor %d  and i have succefully calculate sub_hy",world_rank);    
 }

 
 if((world_rank != 0) && (world_rank != world_size -1) )
 {
   MPI_Recv(&data_recv,1,MPI_FLOAT, world_rank - 1,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  //tag 100 for h_x
   printf("I am  inside processor %d and i have succefully recv data from %d and data is %f \n ",world_rank,world_rank -1,data_recv);
   sub_e_x[0] = sub_e_x[0] + 0.5*(data_recv - sub_h_y[0]);

   data_send = sub_h_y[KE-1];
   MPI_Send(&data_send,1,MPI_FLOAT,world_rank +1 ,100,MPI_COMM_WORLD);    // 100 tag for h_x
   printf("I am  inside processor %d and i have succefully send data to %d and data is %f \n ",world_rank,world_rank +1,data_send);
   for ( k = 1; k < KE; k++ ) {
          sub_e_x[k] = sub_e_x[k] + 0.5*(sub_h_y[k-1] - sub_h_y[k]);
        }

        /* Put a Gaussian pulse in the middle */
        
        if(world_rank == world_size/2)
        {
          T = T+1;
          pulse = exp(-0.5*(pow( (t0-T)/spread, 2.0 )));
          sub_e_x[num_element_per_proc/2] = pulse;
          //printf("xxxxxxxxxxxxxxxx value of  %d and  %d \n", world_rank/2,num_element_per_proc/2);
          //printf("xxxxxxxxxxxxxxxx I am  inside processor %d and writin value %f \n", world_rank,sub_e_x[num_element_per_proc/2]);
          
        }
  
        
	data_send = sub_e_x[0];
//printf("I am  inside processor %d and i have about send data to %d and data is %f \n ",world_rank,world_rank -1,data_send);
        MPI_Send(&data_send,1,MPI_FLOAT,world_rank -1 ,200,MPI_COMM_WORLD);    // 200 tag for e_x
	printf("I am  inside processor %d and i have succefully send data to %d and data is %f \n ",world_rank,world_rank -1,data_send);

//printf("I am  inside processor %d and i have about recv data from %d and data is %f \n ",world_rank,world_rank +1,data_recv);
       // printf("%5.1f \n", sub_e_x[kc]);
        MPI_Recv(&data_recv,1,MPI_FLOAT, world_rank + 1,200,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // 200 tag for e_x
	printf("I am  inside processor %d and i have succefully recv data from %d and data is %f \n ",world_rank,world_rank +1,data_recv);

        /* Calculate the Hy field */
        for ( k = 0; k < KE-1; k++ ) {
          sub_h_y[k] = sub_h_y[k] + .5*(sub_e_x[k] - sub_e_x[k+1] );
        }
        sub_h_y[k] = sub_h_y[k] + 0.5*( sub_e_x[k]- data_recv );

 }




if(world_rank == world_size -1 )
 { 
   MPI_Recv(&data_recv,1,MPI_FLOAT, world_rank - 1,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  //tag 100 for h_x
   printf("I am  inside processor %d and i have succefully recv data from %d and data is %f \n ",world_rank,world_rank -1,data_recv);
   sub_e_x[0] = sub_e_x[0] + 0.5*(data_recv - sub_h_y[0]);

  // data_send = sub_h_y[KE-1];
  // MPI_Send(&data_send,1,MPI_FLOAT,world_rank +1 ,100,MPI_COMM_WORLD);    // 100 tag for h_x
   for ( k = 1; k < KE; k++ ) {
          sub_e_x[k] = sub_e_x[k] + 0.5*(sub_h_y[k-1] - sub_h_y[k]);
        }


        data_send = sub_e_x[0];
        MPI_Send(&data_send,1,MPI_FLOAT,world_rank -1 ,200,MPI_COMM_WORLD);    // 200 tag for e_x
	printf("I am  inside processor %d and i have succefully send data to %d and data is %f \n ",world_rank,world_rank -1,data_send);

       
        /* Calculate the Hy field */
        for ( k = 0; k < KE-1; k++ ) {
          sub_h_y[k] = sub_h_y[k] + .5*(sub_e_x[k] - sub_e_x[k+1] );
        }
        sub_h_y[k] = sub_h_y[k] + 0.5*( sub_e_x[k]- 0.0 );

 }





 MPI_Barrier(MPI_COMM_WORLD); // this is just for test 
 //printf("no of element in processor  = %d \n",  world_rank);

  //float sub_avg = compute_avg(sub_rand_nums,  num_element_per_proc);
  
  // printf("avg in it is %2f \n ",sub_avg);
    for(i =0;i<num_element_per_proc;i++)
   {
     printf("element in processor  = %d is (%2f, %2f) \n ",world_rank,sub_h_y[i],sub_e_x[i]);
   }
 

 if(world_rank == 0)
  {
   final_e = (float*)malloc(num*sizeof(float));
   final_h = (float*)malloc(num*sizeof(float));
  }

 
 MPI_Barrier(MPI_COMM_WORLD);
// MPI_Gather(&sub_avg,1,MPI_FLOAT,sub_avgs,1,MPI_FLOAT,0,MPI_COMM_WORLD);
 
 MPI_Gather(sub_e_x,num_element_per_proc,MPI_FLOAT,final_e,num_element_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Gather(sub_h_y,num_element_per_proc,MPI_FLOAT,final_h,num_element_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
  


  if(world_rank == 0)
  {
   h_y = final_h;
   e_x = final_e;
   //float global_avg =  compute_avg(sub_avgs,  world_size);
   printf("Final Processed value after (%d)th iteration \n",n+1);
   //fprintf( fpEH, "  Start of Time %6d \n", n+1);
    for(i =0;i<num ;i++)
   {
    // fprintf( fpEH, "  %f , %f \n", final_e[i], final_h[i]);
    // fprintf( fpE, "%f,", final_e[i]);
       fprintf( fpH, "%f,", final_h[i]);
     printf(" (%2f, %2f) \n ",final_e[i],final_h[i]);
   }
   //fprintf( fpE, "\n"); 
     fprintf( fpH, "\n"); 
  }

}  // this for outermost for loop time steps

  if(world_rank == 0)
  {
   free(final_e);
   free(final_h);
   //fclose(fpEH);
   //fclose(fpE);
     fclose(fpH);

  }

  free(sub_h_y);
  free(sub_e_x);

  MPI_Barrier(MPI_COMM_WORLD);
 
  MPI_Finalize();




   }

















