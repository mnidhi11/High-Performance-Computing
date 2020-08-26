#include "rk4.h"
#include "pendulum.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//real_t L[10000][3]; // L stores the first periodic point and other periodic points which do not lie in the orbit of any other point in L


//int inner_count; //keeps track of total number of points in L

struct fixedpoint {
   real_t pt[3][2];  // the list of periodic points
   real_t epsilon;   // convergence criterion
   size_t period[3];  // their respective periods
   size_t maxiter;  // maximum number of map iterations
   size_t nfixed;  // maximum number of map iterations
   int p_max;
};



static struct fixedpoint fp = {
   {{-2.407661292643E+00,  6.773267621896E-01},  // period 1
    {-6.099150484926E-01,  1.677463819037E+00},  // period 2
    {-7.218820695570E-01, -4.620944521247E-01}  // period 2
   },  // pt
   1.0E-09,  //  required approach distance to confirm convergence
   {1, 2, 2}, // the respecitve periods
   400, // maxiter
   3, 
   4
};

//  This is one possible way to represent the grid to be used for computing the basin boundary.  

struct grid {
   real_t xmin;  // angular position limits
   real_t xmax;
   real_t ymin;  // angular velocity limits
   real_t ymax;
   size_t resolution;  // number of grid points in each direction
};

// The grid used to generate the basin.  

static struct grid box = {
  MINUSPI, PI,  // xmin, xmax
  MINUSPI, PI,  // ymin, ymax
  40  // resolution
};


// compute_prime_period computes the prime period for every point in the grid 
   



int
compute_prime_period(real_t y[2])
{
     real_t init_value[2]={y[0],y[1]};


      for (int i=0;i<fp.p_max;i++){ //we can iterate for a maximum of p_max times to find the prime period

          poincare(2,y,1);

          if (fabs(init_value[0]-y[0])<(fp.epsilon/2.0) && fabs(init_value[1]-y[1])<(fp.epsilon/2.0)){

  
             return i+1; //returns period =1, 2, 3 or 4

         }

     }
     
     return 0; // return 0 if point is aperiodic

  
}


void do_work(int n,int arr_size,int *basin,int i,real_t L[][3],int row_number,int *inner_count,int me){

   int aperiodic_point=0;

   // points++;

   int count=0;

    real_t temp[arr_size]; // temp stores 2*n elements consisting of position and velocity for the first row in grid 
 

    real_t x=(-PI)+(((PI-(-PI))*(i))/(n-1)); //hardcoded
    
    for (int j_iter=0;j_iter<n;j_iter++){  

        temp[count]=x; //storing position in temp

        count++;

        real_t y=(-PI)+(((PI-(-PI))*(j_iter))/(n-1)); // y denotes velocity

        temp[count]=y; //storing velocity in temp 

        count++;

    }

    poincare(arr_size,temp,fp.maxiter); // step 2 of Algorithm A i.e. Φ j = P (Φ j− 1 ), j = 1 , . . . , M for M=400=fp.maxiter

    for(int t=0;t<arr_size;t+=2){

        real_t temp_array[2]={temp[t],temp[t+1]}; //temp_array stores pair of position and velocity for 2*n elements in the row

        int flag=0; //flag checks if the point lies within the orbit of any point in L

        int result=compute_prime_period(temp_array); // result stores the prime period or 0 if point is aperiodic

        int position=t/2;  // position stores the index of position of an point in the grid

        if(result==0){ // point is aperidoc

           printf("yes\n");

           aperiodic_point++;
             
        }
 
        else
        { 
 
           
                if(*inner_count==0){ 

                    L[*inner_count][0]=temp_array[0];   //storing position

                    L[*inner_count][1]=temp_array[1];   //storing velocity

                    L[*inner_count][2]=result;  //storing prime period  
             
                    *inner_count=(*inner_count)+1;
                    
                 }

                else
                {
                     
                     // step 4 of algorithm

                     for (int iter=0;iter<*inner_count;iter++){

                         int j=0;

                         real_t check_array[2]={L[iter][0],L[iter][1]}; // check_array stores the position and velocity for points in L array which are used to see if the point under consideration lies in the orbit of points in L

                         real_t check_value= check_array[0];

                         while(j<L[iter][2]) { // we iterate Φ M , check distances as before, and repeat up to p max times

                              if(fabs(temp_array[0]-check_value)<fp.epsilon ) { // we check if distance between Φ M and each of the points in L is less than epsilon

                                 *(basin+position+(row_number*n))=iter+1; //storing basin value = 1+ index of point in L 

                                 flag=1;  //distance is less than epsilon.

                                 //Therefore, No need to iterate once more as point lies in the orbit of some point in L

                                 break;

                              }
 
                             poincare(2,check_array,1);

                             check_value=check_array[0];

                             j++;

                         }

                         if(flag==1){ // If one of the distances is less than epsilon ,then mark Φ 0 as being in the basin of the corresponding periodic point in L

                             break;

                         }

                     }

                     if (flag==0){ // Step 5: if distance between Φ M and points in L is not less than epsilon, then store Φ M in L

                         L[*inner_count][0]=temp_array[0];

                         L[*inner_count][1]=temp_array[1];

                         L[*inner_count][2]=result;

                         *inner_count=(*inner_count)+1;
                   
                     }

                 }
         
           
                 

        }


    }
   
  
}

//----------------------------------------------------------------------------
// Compute the basin boundary of each fixed point.  

int*
compute_basin(int n,int nproc,int me,int *root_count,int root_period[],real_t root_period_point[])
{


    int length=n/nproc;

    int row_number=0; 

    int inner_count=0;       

    int arr_size=2*n;

    int *basin = (int *) calloc(((n*n)/nproc), sizeof(*basin));  // initialized to 0

    int i; 

    int i_iter;

    real_t L[length][3];


    if(basin == NULL){

        printf("allocation failure\n");

        return(NULL);

    } 

      // allocation failure
    int a = n/nproc;
    for( i_iter=me*a;i_iter<(me+1)*a;i_iter++) // i_iter is used to iterate over rows. So the rows are distributed amongst k threads)
    {

        i=i_iter;
        do_work(n,arr_size,basin,i,L,row_number,&inner_count,me); // calling do work

        row_number++;
     
    } 

      if(me==0){
        *root_count = inner_count;

        for(int j=0;j<inner_count;j++)
        {
           root_period[j] = (int)L[j][2];
           root_period_point[2*j] = L[j][0];
           root_period_point[2*j+1] = L[j][1];
        }



      }
       
    


      
      if(me!=0){
        int prime_period[inner_count];
        real_t vec_arr[2*inner_count];
        
        for(int k =0;k<inner_count;k++){
            prime_period[k] = (int)L[k][2];          
            vec_arr[2*k]=L[k][0];
            vec_arr[2*k+1]=L[k][1];


       }
           int tag1 = 1;


//MPI_Send(void* data,int count,MPI_Datatype datatype,int destination,int tag,MPI_Comm communicator)//
           MPI_Send(&inner_count,1,MPI_INT,0,tag1,MPI_COMM_WORLD);                        //Sends the number of fixed points found by non root process to root process

           int tag2 = 2;

           MPI_Send(prime_period,inner_count,MPI_INT,0,tag2,MPI_COMM_WORLD);              //Sends the prime periods of fixed points found by non root process to root process                       

           int tag3=3;

           MPI_Send(vec_arr,2*inner_count,MPI_DOUBLE,0,tag3,MPI_COMM_WORLD);              //Sends the fixed points found by non root process to root process


           int tag4=4;

           MPI_Send(basin,(n*n)/nproc,MPI_INT,0,tag4,MPI_COMM_WORLD);                     //Sends the basin found by non root process to root process

          }

     
    return(basin);

}
  



char
 check_point(real_t y[])
 {
  real_t pB[] ={-6.099150484926E-01,  1.677463819037E+00};
  real_t pC[] ={-7.218820695570E-01, -4.620944521247E-01};

  poincare(2,pB,1);
  poincare(2,pC,1);

     


    if((fabs(y[0]-fp.pt[0][0])<=fp.epsilon) && (fabs(y[1]-fp.pt[0][1])<=fp.epsilon)) //Compute distance between (phi(M)) and fixed point A
     {
    return 'A';
    }
    
      else if((fabs(y[0]-fp.pt[1][0])<=fp.epsilon)  && (fabs(y[1]-fp.pt[1][1])<=fp.epsilon)) //Compute distance between (phi(M)) and fixed point B
    {
      
      return 'B';
    }
    else if((fabs(y[0]-pB[0])<=fp.epsilon) && (fabs(y[1]-pB[1])<=fp.epsilon))  //Compute distance between (phi(M)) and point P(B)
    {
    
      return 'B';
    }
  
    else if((fabs(y[0]-fp.pt[2][0])<=fp.epsilon)  && (fabs(y[1]-fp.pt[2][1])<=fp.epsilon))  //Compute distance between (phi(M)) and fixed point C
     { 
     
      return 'C'; 
    }

    else if((fabs(y[0]-pC[0])<=fp.epsilon) && (fabs(y[1]-pC[1])<=fp.epsilon)) //Compute distance between (phi(M)) and point P(C)
    {
      return 'C';
      
    }
else{
  exit(0);
}
    
  }



void modify_grid(real_t root_period_point[],int grid_size,int period[],int L_count,int basin_array[],int proc_size,int sender,int total_basin_array[],int grid_iter,int n,real_t period_point[])
{
  
  //Logic for post processing 

int arr[L_count+1];
char root_arr[L_count];
char period_val[L_count];
real_t y[2];
real_t r[2];  


for(int iter=0;iter<L_count;iter++)                                      //Stores the fixed points found by root process to temporary array
{
r[0] = root_period_point[2*iter];
r[1] = root_period_point[2*iter+1];                                                   
char temp_r = check_point(r);
root_arr[iter] = temp_r;
}


  for(int i = 0;i<L_count;i++)                                          //Loop for non root processes
  {
  y[0] = period_point[2*i];                                             //Stores values of fixed points found by non root processes to temporary array         
  y[1] = period_point[2*i+1];

    period_val[i] = check_point(y);                                                                                    

   for(int j = 0;j<L_count;j++)                                         //Loop for root process
   {
   if(root_arr[j] == period_val[i])
    
    {
 
      arr[j+1] = i+1;                                                   //Reference array to modify the basin values
      
      break;
    }
   }
   
  }

  for(int iter=0;iter<proc_size;iter++)
  {
    for(int iter_i=1;iter_i<=L_count;iter_i++)
    {
      if(basin_array[iter] == arr[iter_i])                              //Modify the basin values according to the reference array.
      {
        basin_array[iter] = iter_i;
      
        break;

    }
     
  }
 
}

}


//----------------------------------------------------------------------------
// Main program.  
// The version below simply allocates an NxN array of integers for the
// basin boundary array.

int
main(int argc, char **argv)
{
  

    int me, nproc,n;
    int L_count;
    int sum = 0;
    int index=0;
    int index_point = 0;
    int basin_index = 0;

    char const *file = "input.dat";
  


    MPI_Status status;
    real_t param[3];

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    if(me == 0) {

      FILE *inputfp = fopen(file,"r");

      if(inputfp==NULL){
        printf("Cannot open file\n");
      }

      for(int i=0;i<3;i++)
      {
        fscanf(inputfp,"%lf",&param[i]);
      }


    }

    MPI_Bcast(param, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);                                     //Broadcast delta,Force and number of processes to non root processes

    delta = param[0];

    force = param[1];

    n =  (int)param[2];

    int root_count=0;
    int root_period[n*n];
    real_t root_period_point[2*n*n];
    
    int sender;
    int total_basin_array[n*n];
    int *basin;
    int grid_size = n*n;
    int proc_size = (n*n)/nproc;
    int grid_iter=n/nproc;
    real_t period_point_orbit[100][3];
    

    if(me==0)
    {
   
                                                                                                            
    basin = compute_basin(n,nproc,me,&root_count,root_period,root_period_point);              //Calculate basin for root process
    sum = sum + root_count;
    for(int j=0;j< root_count;j++)
    {
         
     period_point_orbit[index][2] = root_period[j]; 
     index++;

    period_point_orbit[index_point][0] = root_period_point[2*j];
    period_point_orbit[index_point][1] = root_period_point[2*j+1];

    index_point++;

    }
    
    for(int j = 0;j<(n*n)/nproc;j++)                                                      //Store basin values found by root process
     {
     total_basin_array[basin_index] = basin[j];
     basin_index++;



      }
    

    }

      if(me!=0)                                                                           //Calculate basin of non root processes
      {
        basin = compute_basin(n,nproc,me,&root_count,root_period,root_period_point);

      }
   
 
    
    if(me==0){
   
    for(int p=0;p<nproc-1;p++)
    {
      MPI_Recv(&L_count,1,MPI_INT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status); //Receives total count of periods

      sender = status.MPI_SOURCE;
    
      sum= sum + L_count;

     
      int period[L_count];

      MPI_Recv(period,L_count,MPI_INT,sender,2,MPI_COMM_WORLD,&status);   //Receives period point array 
   

      for(int loop =0;loop<L_count;loop++)
      { 
        period_point_orbit[index][2] =  period[loop];                                                    //Stores period of fixed points
        index++;
      }
  
      real_t period_point[2*n];
   
      MPI_Recv(period_point,2*L_count,MPI_DOUBLE,sender,3,MPI_COMM_WORLD,&status);  //Receives an array with periodic point 
   

    
        for(int i=0;i<L_count;i++)
        {
        
        period_point_orbit[index_point][0] =  period_point[2*i];                    //Stores periodic point
        period_point_orbit[index_point][1] =  period_point[2*i+1]; 
        index_point++;
       }
   

       int basin_array[(n*n)/nproc];
    
       MPI_Recv(basin_array,(n*n)/nproc,MPI_INT,sender,4,MPI_COMM_WORLD,&status);  //Receives basin array
  

      modify_grid(root_period_point,grid_size,period,L_count,basin_array,proc_size,sender,total_basin_array,grid_iter,n,period_point);

   
    int count=0;
    for(int i_iter=sender*(grid_iter);i_iter<(sender+1)*(grid_iter);i_iter++)
     {
       for(int j_iter=0;j_iter<n;j_iter++)
       {
          total_basin_array[j_iter+(i_iter*n)]=basin_array[count];     //combines basin array of all processes
          count++;
        }
     }
    }
  
  }


   


      MPI_Barrier(MPI_COMM_WORLD);

    // Write out the results in a binary format compatible with the MATLAB
    // script wada.m for visualization.
      if(me==0){
        
 
  int iter_j;  

for(int iter_i=0; iter_i<index_point; iter_i++)
    {
        for(iter_j=iter_i+1; iter_j<index_point; iter_j++)
        {
            /* If any duplicate found */
            if(fabs((period_point_orbit[iter_i][0]-period_point_orbit[iter_j][0])<=fp.epsilon) && fabs((period_point_orbit[iter_i][1]-period_point_orbit[iter_j][1])<=fp.epsilon))
            {
                
                /* Delete the current duplicate element */
                for(int iter_k=iter_j; iter_k<index_point; iter_k++)
                {   
                    period_point_orbit[iter_k][0]=period_point_orbit[iter_k+1][0];
                    period_point_orbit[iter_k][1]=period_point_orbit[iter_k+1][1];
                    period_point_orbit[iter_k][2]=period_point_orbit[iter_k+1][2];
                }
                /* Decrement size after removing duplicate element */
                index_point--;
                /* If shifting of elements occur then don't increment j */
                iter_j--;
            }
        }

       

    }


  
        for(int i=0;i<index_point;i++)
      {
        printf("Periodic points are: %g ,%g with periods %g\n",period_point_orbit[i][0],period_point_orbit[i][1],period_point_orbit[i][2]);
      }

    const char *filename = "basin.dat";

    FILE *outfp = fopen(filename, "wb");

    if(outfp == NULL) {  // cannot open output file, so halt

       perror(filename);

       return(EXIT_FAILURE);

    }

    fwrite(&box, sizeof(box.xmin), 4, outfp);  // xmin through ymax

    fwrite(&fp.epsilon, sizeof(fp.epsilon), 1, outfp);

    fwrite(&n, sizeof(n), 1, outfp);

    fwrite(total_basin_array,sizeof(int) ,n*n, outfp);  // N x N basin boundary

    fclose(outfp);

  //  free(basin);
  }

    MPI_Finalize();

    //return(EXIT_SUCCESS);
       return 0;
}
