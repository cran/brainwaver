/* Study of the connexity of networks */


/* Internal function to return the component reachable by the 
 node  num_node*/

void reachable_component(int num_node, int nnod, int *conn, int *rc, double *charac_path_length)
{

  int v,j,minl;
  int count=0;
  int l[nnod];
  int tmp=0;

    v = num_node;  /* analysed node is initially num_node */

    /* initialise rc and l */
    
    for(j = 0; j < nnod; j++)  
    {
      if(j == v)  /* distance and inclusion of initial vertex in tree */
      {
        rc[j] = 1; 
        l[j] = 0;
      }
      else  /* initially we don't now shortest path to any vertex */
      {
        rc[j] = 0;
        l[j] = 1e5;
      }  
    }
    

    
    do
    {
      /* update distances (labels) of arcs containing v and a node not included in rc */
    
      for(j = 0; j < nnod; j++)
      {
        if((rc[j] == 0) && (conn[ij(v,j,nnod)] == 1))  /* j is node not in rc and has arc with v */
        {
          if(l[v] + 1 < l[j])  /* we should update label of j */
          {
            l[j] = l[v] + 1;
          }
        }
      }
     
      /* look for new v */
    
      minl = 1e5;
      for(j = 0; j < nnod; j++)  
      {
        if((rc[j] == 0) && (l[j] < minl))
        {
          minl = l[j];
          v = j;
        }
      }
    
      if(minl < 1e5)  /* there is a new vertex to include in T? */
        {
	  rc[v] = 1;
	  tmp+=l[v];
	  count++;
	}

    }    
    while(minl < 1e5);  /* while there is a new vertex v included in T */  

    if(count!=0) 
      *charac_path_length = (double) (tmp) / (double) (count);
    else 
      *charac_path_length=-1.0;
    
}

/*R functions*/


void Rreachable_component(int *pnum_node, int *pnnod, int *conn, double *Rrc)
{
  int j;
  double charac_path_length;
  int num_node=*pnum_node-1;
   int nnod=*pnnod;
   int rc[nnod];
    reachable_component(num_node,nnod,conn,rc,&charac_path_length);    
    //printf("charac=%e \n",charac_path_length);
    for(j = 0; j < nnod; j++) Rrc[j]=(double) rc[j];
    Rrc[nnod]=charac_path_length;

}


void Rconnex_components (int *pnnod, int *conn, double *acc)
{
  int j,i,max;
  int nnod=*pnnod;
  int i_max=-1;
  int cc_max=0;
  int num_cc=0;
  int rc[nnod];
  double l_tmp;

  /* Initialisation of acc */

  for(j = 0; j < nnod; j++)
    { 
      acc[j]=-1;
    }

  i=0;

  do
    {
      if(acc[i]==-1)
	{
	  reachable_component(i,nnod,conn,rc,&l_tmp);
	  max = 0;
	  for(j = 0; j < nnod; j++)
	    {
	      if(rc[j]==1) 
		{
		  acc[j]=i+1;
		  max++;
		}
	    }
	  if(cc_max<max)
	    {
	      cc_max=max;
	      i_max=i+1;
	    }
	  num_cc++;	
	}
      i++;
    
    }while(i<nnod); 


  acc[nnod]=i_max;
  acc[nnod+1]=cc_max;
  acc[nnod+2]=num_cc;
}


void Rcharac_path_length(int *pnum_node, int *pnnod, int *pcc_size, int *conn, double *Lp)
{
  int j;
  int nnod=*pnnod;
  int num_node=*pnum_node-1;
  int cc_size=*pcc_size;
  int rc[nnod];
  int rc_tmp[nnod];
  double l;

  reachable_component(num_node,nnod,conn,rc,&l);
  *Lp=l;
  //printf("num %i, carac path length = %e \n",num_node,l); 

    for(j=0; j < nnod; j++)
      {
	if((rc[j]==1)&&(j!=num_node)){
	  reachable_component(j,nnod,conn,rc_tmp,&l);
	  // printf("num %i, carac path length = %e \n",j,l); 
	  *Lp+=l;
	}
      }

    *Lp/=(double) (cc_size);

 
  
}

      









