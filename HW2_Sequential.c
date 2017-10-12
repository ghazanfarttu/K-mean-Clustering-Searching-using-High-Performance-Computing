/*
 * =====================================================================================
 *   Homework 2:  K-Means Clustering and Searching using HPC
 *   Filename:  HW2.c
 *
 *   Author:  Ghazanfar Ali, ghazanfar.ali@ttu.edu
 			  Rishika,							 
 * =====================================================================================
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>

double distance(int dims, double *data, double *cluster)
{
  int i; double diff=0.0, dist = 0.0;
  for (i=0; i<dims; i++)
    {
      diff = data[i] - cluster[i];
      dist += diff * diff;
    }
  return dist;
}

int assignment_change_count(int n, int *a, int *b)
{
    int change_count = 0;

    for (int ii = 0; ii < n; ii++)
      if (a[ii] != b[ii])
        change_count++;
        
    return change_count;
}

int calculate_K_radius(double *data, int dim, int proc_data, int K, int *cluster_assign, double *centroids, double *cluster_radius)
{
	int d, i, j, n, x, max_index;
	double *member_dist, max_dist = 0.0;
	double tempData[dim], singleCentroid[dim];
	
	member_dist = (double*) malloc(proc_data * sizeof(double));
    
	
	// Calculate the  K radius
	for (i = 0; i < K; i++ )
	{
		
		for (d = 0; d < proc_data; d++)
			member_dist[d] = 0.0;
		
		for (n =0; n < dim; n++)
			singleCentroid[n] = centroids[i*dim + n];
				
		for (d = 0; d < proc_data; d++)
		{
	  		if (cluster_assign[d] != i)					
	  			continue;
	
	  		for (n = 0; n <dim; n++)
	    		tempData[n] = data[d*dim + n];
     
	  		// find distance between clusters and its members
	  		member_dist[d] = distance(dim, tempData, singleCentroid);
	  
		}//data loop
		
		// find the max of min distance array
		max_dist = member_dist[0];
	 		
	 	for (x=1; x < proc_data; x++)
	 	{
			if (cluster_assign[x] != i)					
	  			continue;
	  							
			if ( max_dist < member_dist[x])
	 		{
	 			max_dist = member_dist[x];
	 			max_index = x;
	 		}
	 		
	 	}
	 	// copy the furthest data point
	 	for (n = 0; n <dim; n++)
	    	tempData[n] = data[max_index*dim + n];
	    // calculate distance between max data point and corresponding cluster
	 	cluster_radius[i] = distance(dim, tempData, singleCentroid);
	 	
	 	/*	
		for (n = 0; n < dim; n++)
			cluster_radius[(i*dim + n)] = abs ((data[(max_index * dim + n)] - centroids[(i*dim + n)]));	
		*/
	} //K-Radii		
}

int kmean_search(double *query_point, double *data, int dim, int proc_data, int K, int *cluster_assign, double *centroids, int *cluster_size, double *cluster_radius, double *nearest_point)
{
  	// find nearest cluster
  	double max_dist, dist, refDist, singleCentroid[dim], min_dist[proc_data], tempData[dim], d_clusters[K];

  	int k, i, j, n, x, d;
	
  	refDist = DBL_MAX;
	
	for (d = 0; d < proc_data; d++)
		min_dist[d] = 0.0;
		
  	for (j =0; j < K; j++)
  	{
	  
		for (n =0; n < dim; n++)
	  		singleCentroid[n] = centroids[j*dim + n];
	      
		// find distance
		dist = distance(dim, query_point, singleCentroid);
		//printf("large \n");
	  
		if ( dist < refDist )
 		{
			k = j;
			refDist = dist;
			// printf("less \n");
		}
		// distance of query point to j-th cluster
		d_clusters[j] = dist;
  	}// cluster loop
  	 	
	//printf ("\n nearest cluster - k = %d \n", k);
	//printf ("\n dist =  %.0lf \n",refDist);
	
	for (n =0; n < dim; n++)
		singleCentroid[n] = centroids[k*dim + n];
	
	// find nearest point in the cluster
	for (d = 0; d < proc_data; d++)
	{
	  if (cluster_assign[d] != k)					
	  	continue;
	
	  for (n = 0; n <dim; n++)
	    tempData[n] = data[d*dim + n];
     
	  // find distance between clusters and its members
	  min_dist[d] = distance(dim, tempData, singleCentroid);
	  		
	}//data loop
	/*
	for (n=0; n < K; n++)
		printf (" d_clusters[%d] =  %.0lf   ",n,d_clusters[n]);
	
	printf ("\n");
	*/
	/*
	for(i=0; i < proc_data; i++)
	{
		
		if (min_dist[i] == refDist)
		{
			x = i;
			break;
		}
	}
	for (i = 0;i < dim; i++)
		nearest_point [i] = data[x * dim + i];
	return 0;
	*/
	
	// find the nearest point to the query point
	double smallestDiff = refDist - min_dist[0];
	int closest = 0; //index of the current closest number
	double currentDiff;
	
	for (i = 1; i < proc_data; i++) 
	{
		if (cluster_assign[i] != k)					
	  		continue;
	  		
  		currentDiff = refDist - min_dist[i];
  		
  		if (currentDiff == 0) 
  		{
  			// query point is in the data set
    		closest = i;
  			break;
  		}
  		else if (currentDiff <= smallestDiff) 
  		{
  			// query point is not in the data set
    		smallestDiff = currentDiff;
    		closest = i;
  		}
	}
	//printf (" NEAREST POINT = %d   ", closest);
	// copy the nearest point
	for (i = 0;i < dim; i++)
		nearest_point [i] = data[closest * dim + i];
		
	// check the plausible closest points in other clusters
	int closest_pt_ind[K];
	double qp_np_dist [K],qp_dist;
	
	for (i = 0; i < K; i++)
	{
		closest_pt_ind[i] = -1;
		qp_np_dist [i] = -99999.99;
	}
	
	closest_pt_ind[k] = closest;
	dist = distance (dim, nearest_point, query_point);
	qp_np_dist[k] = dist;
	  	
	for (x = 0; x < K; x++)
	{	
		if ( (d_clusters[x] - cluster_radius[x]) <  dist)
			{
				// find the closest point in cluster x
				//printf ("\n CLOSET POINT IN OTHER CLUSTERS POSSIBLE \n");
							
				// extract the x centroid
				for (n =0; n < dim; n++)
					singleCentroid[n] = centroids[x*dim + n];
		
				// find nearest point in the cluster
				for (d = 0; d < proc_data; d++)
				{
	  				if (cluster_assign[d] != x)					
	  					continue;
	
	  				for (n = 0; n <dim; n++)
	    				tempData[n] = data[d*dim + n];
     
	  				// find distance between cluster x and all its members
	  				min_dist[d] = distance(dim, tempData, singleCentroid);
	  		
				}//data loop
	
				// find the nearest point to the query point
				smallestDiff = d_clusters[x] - min_dist[0];
				closest = 0; //index of the current closest number
			
				for (i = 1; i < proc_data; i++) 
				{
					if (cluster_assign[i] != x)					
	  					continue;
	  		
  					currentDiff = d_clusters[x] - min_dist[i];
  					if (currentDiff == 0) 
  					{
  						// query point is in the data set
    					closest = i;
  						break;
  					}
  					else if (currentDiff < smallestDiff) 
  					{
    					smallestDiff = currentDiff;
    					closest = i;
  					}
				}
	
				// copy the nearest point
				for (i = 0;i < dim; i++)
					nearest_point [i] = data[closest * dim + i];
			
				closest_pt_ind[x] = closest;
				qp_np_dist[x] = distance (dim, nearest_point, query_point);		
			} // if 
	} // clusters
	
	// find the "real" closest point
		int min_index; double min_dis;
		min_dis = qp_np_dist[0];
	 	min_index = 0;	
	 	for (x=1; x < K; x++)
	 	{
	 		//printf ("\n  qp_np_dist[%d]=%f\n", x-1, qp_np_dist[x-1]);	
	 		if (qp_np_dist[x] == -99999.99)
	 			continue;		
			if (qp_np_dist[x] == 0)
			{
				min_index = x;
				break;
			}
			if ( min_dis > qp_np_dist[x])
	 		{
	 			min_dis = qp_np_dist[x];
	 			min_index = x;
	 		}
	 		
	 	}
	 	//printf ("\n final cluster = %d\n", min_index);
	 	
	 	min_index = closest_pt_ind[min_index];
	 	// copy the furthest data point
	 	//printf ("\n nearest point is = %d\n", min_index);
	 	for (n = 0; n <dim; n++)
	    	nearest_point[n] = data[min_index*dim + n];
}

int build_K_centroids(int dim, int proc_data, double *data, int K, double *centroids)
{
	int d, i, j, n, x,y=0, max_index, already_taken[K], t, flag;
	double *min_dist, refDist=0.0, dist = 0.0, max_dist = 0.0;
	double tempData[dim], singleCentroid[dim];
	
	min_dist = (double*) malloc(proc_data * sizeof(double));
	
	//get the first K Clusters elements from data[] as the initial cluster centers
    srand(time(NULL));
	int random_index = rand()%10000;
	//printf("Index = %d\n", random_index);
    for (i=0; i < dim; i++)
	centroids[i] = data[random_index * dim + i];
   	already_taken[y++] = random_index;
	
	// Construct K-1 centroids
	for (i = 1; i < K; i++ )
	{
				
		flag = 0;
		for (d = 0; d < proc_data; d++)
		{
						
	  		for (n = 0; n <dim; n++)
	    		tempData[n] = data[d*dim + n];
     
	  		refDist = DBL_MAX;
	
	  		for (j =0; j < i; j++)
	    	{
	  
	      		for (n =0; n < dim; n++)
					singleCentroid[n] = centroids[j*dim + n];
	      
	      		// find distance
	      		dist = distance(dim, tempData, singleCentroid);
	  
	      		if ( dist < refDist )
 				{
		  			//k = j;
		  			refDist = dist;
		  
				}
	     
	    	} // cluster loop
			// copy the mini distance of data [d] in min distance array
			min_dist[d] = dist;
			
		}//data loop
		
		// find the max of min distance array
		max_dist = min_dist[0];
	 		
	 		for (x=1; x < proc_data; x++)
	 		{
	 			flag = 0;
	 			for (t=0; t<y; t++)
					if (already_taken[t] == x)
						flag = 1;
						
							
				if ( (max_dist < min_dist[x]) && flag == 0)
	 			{
	 				max_dist = min_dist[x];
	 				max_index = x;
	 			}
	 		
	 		}
	 		
	 	already_taken[y++] = max_index;	
		for (n = 0; n < dim; n++)
			centroids[i*dim + n] = data[max_index * dim + n];
	} // k cluster construction loop
	
}

int chkHandleEmptyCluster(int K,int dim, int proc_data,  double *data, double *tmpCentroids,  int *tmpClusterSize, int *clusterAssign)                                
{
  
  int k, i, j, f, s, m, prev_cluster=0, flag=0;
  double dist =0.0, refDist, tempData[dim], singleCentroid[dim];

  int maxDataInd;
  
  int d,  n, x,y=0, already_taken[K], t;
	double *min_dist, max_dist = 0.0;
	
	
	min_dist = (double*) malloc(proc_data * sizeof(double));

  for (k=0; k < K; k++)
  {
    
    if ( tmpClusterSize[k] == 0 )
	{
	  flag = 0;
	  for (d = 0; d < proc_data; d++)
		{
						
	  		for (n = 0; n <dim; n++)
	    		tempData[n] = data[d*dim + n];
     
	  		refDist = DBL_MAX;
	
	  		for (j =0; j < K; j++)
	    	{
	  
	      		for (n =0; n < dim; n++)
					singleCentroid[n] = tmpCentroids[j*dim + n];
	      
	      		// find distance
	      		dist = distance(dim, tempData, singleCentroid);
	  
	      		if ( dist < refDist )
 				{
		  			//k = j;
		  			refDist = dist;
		  
				}
	     
	    	} // cluster loop
			// copy the mini distance of data [d] in min distance array
			min_dist[d] = dist;
			
		}//data loop
		
		// find the max of min distance array
		max_dist = min_dist[0];
	 		
	 		for (x=1; x < proc_data; x++)
	 		{
	 			flag = 0;
	 			for (t=0; t<y; t++)
					if (already_taken[t] == x)
						flag = 1;
						
							
				if ( (max_dist < min_dist[x]) && flag == 0)
	 			{
	 				max_dist = min_dist[x];
	 				maxDataInd = x;
	 			}
	 		
	 		}
	 		
	 	already_taken[y++] = maxDataInd;	
		for (n = 0; n < dim; n++)
			tmpCentroids[i*dim + n] = data[maxDataInd * dim + n];
	  // set empty cluster id to  'cluster_assign' and increase the size                                                                                                         
	  prev_cluster = clusterAssign [maxDataInd];
	  //clusterAssign[i] = maxDataInd;
	  clusterAssign [maxDataInd] = k;
	  tmpClusterSize[k]++;
	  
	  
	  // decrease the size of previous cluster by 1
	  tmpClusterSize[prev_cluster]--;

	  // assign furthest data point as centroid of the empty cluster                                                                                                                      
	  for (m=0; m<dim; m++)
	  	tmpCentroids[(k*dim)+m] += data[(maxDataInd*dim) + m];
	
	 // recalculate the centroid of the previous cluster after removing furthest point
	 for (m=0; m<dim; m++)
	 {
	 	tmpCentroids[(prev_cluster*dim)+m] *= tmpClusterSize[prev_cluster];
		tmpCentroids[(prev_cluster*dim)+m] -= data[(maxDataInd*dim) + m];
    	tmpCentroids[(prev_cluster*dim)+m] /= --tmpClusterSize[prev_cluster];
	 }  
	}// if					\
      
  }// clusters for loop			\
  
} // external
       
int kmean(double *data, int dim, int proc_data, int K, int *cluster_assign, double *centroids, int *cluster_size, double *cluster_radius)
{
  double dist, refDist, tempData[dim], singleCentroid[dim];
 
  double *tmpCentroids;
  tmpCentroids = (double*) malloc((K * dim) * sizeof(double));

  int *tmpClusterSize, *prev_cluster_assign;
  tmpClusterSize = (int*) malloc(K * sizeof(int));
  prev_cluster_assign = (int*) malloc(proc_data * sizeof(int));

  int k, i, j, m, n,s, iterate, y, flag = 0, change_count;

  //copy initial clusters values into tmpClusters
  for(i=0;i<K*dim; i++)
    tmpCentroids[i] = centroids[i];
  
  while (1)
    {
			
	// copy the initial cluster membership and set to zero
      for (i =0; i < proc_data; i++)
	  {	
	  	prev_cluster_assign [i] = cluster_assign[i];
		cluster_assign[i] = 0; 
 	  }
 	  
      for (i =0; i < K; i++)
		tmpClusterSize[i] = 0;

	

      for (i = 0; i < proc_data; i++)
	{
	  for (n = 0; n <dim; n++)
	    tempData[n] = data[i*dim + n];
     
	  refDist = DBL_MAX;
	
	  for (j =0; j < K; j++)
	    {
	  
	      for (n =0; n < dim; n++)
			singleCentroid[n] = tmpCentroids[j*dim + n];
	      
	      // find distance
	      dist = distance(dim, tempData, singleCentroid);
	      //printf("large \n");
	  
	      if ( dist < refDist )
 		{
		  k = j;
		  refDist = dist;
		  // printf("less \n");
		}
	     
	    } // cluster
	
	  tmpClusterSize[k]++;
          cluster_assign[i] =k;

	}//data loop


      // calculate mean of all datapoints assigned to corresponding clusters 
      int a, b;
      for (a = 0; a < K; a++)
	for (b=0; b < dim; b++)
	  tmpCentroids[a*dim + b] = 0.0;
  
      int  u, w,dt, size;

      for (s =0; s < K; s++)// for all clusters K
	{
	  size = tmpClusterSize[s];// find members for a cluser k
      
	  if (size>0)
	    {

	      for (u = 0; u < proc_data; u++)
		{
		  if (cluster_assign[u] == s)
		    {
		      for (w=0;w<dim;w++)
			tmpCentroids[s*dim+w] += data[u*dim+w]; 
		    } //end of member
		} //end of data
	      for (dt=0;dt<dim;dt++) // divide each centroid by number of data points contained in the cluster
		tmpCentroids[s*dim + dt] /= size;
	    } // end of outer if
	} // end of outer for                      

      // Check and handle empty clusters
      chkHandleEmptyCluster(K,dim, proc_data,  data,  tmpCentroids, tmpClusterSize, cluster_assign);

	  change_count = assignment_change_count(proc_data, cluster_assign, prev_cluster_assign);
	  printf ("change count = %d\n", change_count);
	  if ( change_count == 0)
		break;
    } // while loop
    
    // replace initial centroids with final values
    for (i=0;i<K*dim;i++)
    	centroids [i] = tmpCentroids[i];
  		
  calculate_K_radius(data, dim, proc_data, K, cluster_assign, centroids, cluster_radius);
  /*
  printf("\n ...RADII...\n");	
  
  for (i = 0; i < K; i++)
  	printf("   r[%d] = %.0lf   ", i, cluster_radius[i]);
  	
  printf("\n ...RADII...\n");	
  */
  free (tmpClusterSize);
  free (tmpCentroids);
}

int main(int argc, char** argv) 
{

int rank, procSize;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &procSize);
  
  // ELAPSED TIME
  double secs = 0.0, totalClockTime = 0.0;
  struct timeval start, stop;
  gettimeofday(&start, NULL);

  // CPU Clock TIME
  clock_t begin, finish;
  begin = clock();

  //printf ("\nNumber of procs = %d\n", procSize);
  
  //10000/128
  int ndata = 10000, dim = 128, i,j=0, n, K, seeds[4]={1,2,3,4};
   int seed_index, proc_data;
   
   //K = pow(2,ceil (log(sqrt(ndata))));
	K = ceil (sqrt(ndata));
  proc_data = ndata/procSize;
                                                                                                                       
  int *cluster_assign;
  cluster_assign = (int*) malloc(proc_data * dim * sizeof(int));                                                                                                                     
                                                                                                                                              
  double *centroids;
  centroids = (double*) malloc((K * dim) * sizeof(double)); 
                                                                                                                   
  double *data;
  data  = (double*) malloc(proc_data * dim * sizeof(double));
                                                                                                                                                                                
  int *cluster_size;
  cluster_size = (int*) malloc(K * sizeof(int)); 
  
  double *cluster_radius;
  cluster_radius = (double*) malloc((K * dim) * sizeof(double));


//every process generates its own data
  
// rank starting index in the seed array  
  seed_index=rank * 4/procSize;
                                                                                                         
  // seed = seeds[rank*numSeeds]; // calculate starting point of seeds[]  for each rank
  for (i=0; i < proc_data * dim; i++)
    {
      //320000
      if ((i)%320000 == 0)
	  	srand (seeds[seed_index++]);
	
	
      data[i] = (rand());
      data[i] = (((data[i] / RAND_MAX) * 100.0) + 1);
      /*printf (" %.0lf ",data[i]);
     	if (i%50==0)
     		printf("\n"); */
    }
   //printf("\n"); 
   
	printf ("\nK = %d\n", K);
  
  // constructing k-means build_K_means
  build_K_centroids(dim, proc_data, data, K, centroids);
 
  kmean(data, dim, proc_data, K, cluster_assign, centroids, cluster_size, cluster_radius);
  
  double *nearest_point/*, *query_point*/;
  nearest_point  = (double*) malloc(dim * sizeof(double));
  //query_point = (double*) malloc(dim * sizeof(double));
 
  double query_point[128] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
          18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
          35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
          52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
          69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85,
          86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 10, 11, 22,
          33, 44, 55, 66, 77, 88, 99, 10, 11, 12, 13, 14, 15, 16,
          17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
	/*	
  //get a random number from the data set
	
	srand(time(NULL));
	int random_index = rand()%10000;
	//printf("Point Index = %d\n", random_index);
    
    for (i=0; i < dim; i++)
		query_point[i] = data[random_index * dim + i];
	*/
  kmean_search(query_point, data, dim, proc_data, K, cluster_assign, centroids, cluster_size, cluster_radius, nearest_point);
  
  printf("\n*********** Q U E R Y          P O I N T ***********\n");
  
  for (i = 0; i < dim; i++)
	printf (" %.0lf ",query_point[i]);
		
  printf("\n");
  	
  printf("\n*********** N E A R E ST        P O I N T ***********\n");
  
  for (i = 0; i < dim; i++)
	printf (" %.0lf ",nearest_point[i]);
		
  printf ("\n");
  
	
 free (cluster_size);
  free (centroids);
  // close time
  
  finish = clock();

  totalClockTime = ((double)(finish - begin)/CLOCKS_PER_SEC);
  printf ("\nCPU clock time = %f\n", totalClockTime);

  gettimeofday(&stop, NULL);
  secs = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
  printf ("\nElapsed Time = %f\n",secs);

	MPI_Finalize();
	
  return 0;
}
