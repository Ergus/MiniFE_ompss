//@HEADER
// ************************************************************************
//
// MiniFE: Simple Finite Element Assembly and Solve
// Copyright (2006-2013) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// ************************************************************************
//@HEADER

#include <stdio.h>
#include <stdlib.h>

#include <Box.hpp>
#include <BoxPartition.hpp>

/*--------------------------------------------------------------------*/

static int box_map_local_entry( const Box& box ,
                                const int ghost ,
                                int local_x ,
                                int local_y ,
                                int local_z )
{
  const int nx = 2 * ghost + box[0][1] - box[0][0] ;
  const int ny = 2 * ghost + box[1][1] - box[1][0] ;
  const int nz = 2 * ghost + box[2][1] - box[2][0] ;
  int result = -1 ;

  local_x += ghost ;
  local_y += ghost ;
  local_z += ghost ;

  if ( 0 <= local_x && local_x < nx &&
       0 <= local_y && local_y < ny &&
       0 <= local_z && local_z < nz ) {

    result = local_z * ny * nx + local_y * nx + local_x ;
  }
  return result ;
}

int box_map_local( const Box& box_local,
                   const int ghost ,
                   const int box_local_map[] ,
                   const int local_x ,
                   const int local_y ,
                   const int local_z )
{
  int result = box_map_local_entry(box_local,ghost,local_x,local_y,local_z);

  if ( 0 <= result ) {
    result = box_local_map[ result ];
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/* Recursively split a box into into (up-ip) sub-boxes */

void box_partition( int ip, int up, int axis, const Box &box, Box *p_box )
{
	const int np = up - ip ;

	if ( np == 1 ) {
		assert(box.get_num_ids() != 0);
		p_box[ip] = box;
		return;
	}

	const int n = box[axis][1] - box[axis][0] ;
	const int np_low = np / 2 ;  /* Rounded down */
	const int np_upp = np - np_low ;

	const int n_upp = (int) (((double) n) * ( ((double)np_upp) / ((double)np)));
	const int n_low = n - n_upp ;
	const int next_axis = ( axis + 2 ) % 3 ;

	if ( np_low ) { /* P = [ip,ip+np_low) */
		Box dbox(box) ;

		dbox[axis][1] = dbox[ axis ][0] + n_low ;

		box_partition(ip, ip + np_low, next_axis, dbox, p_box);
	}

	if ( np_upp ) { /* P = [ip+np_low,ip+np_low+np_upp) */
		Box dbox(box);

		ip += np_low ;
		dbox[axis][0] += n_low ;
		dbox[axis][1]  = dbox[axis][0] + n_upp ;

		box_partition(ip, ip + np_upp, next_axis, dbox, p_box);
	}
}

/*--------------------------------------------------------------------*/

static int box_disjoint(const Box &a , const Box &b)
{
  return a[0][1] <= b[0][0] || b[0][1] <= a[0][0] ||
         a[1][1] <= b[1][0] || b[1][1] <= a[1][0] ||
         a[2][1] <= b[2][0] || b[2][1] <= a[2][0] ;
}

static void resize_int(int **a , int *allocLen , int newLen)
{
	int k = 32;
	while ( k < newLen )
		k <<= 1 ;

	if ( *a == NULL )
		*a = (int*)malloc( sizeof(int)*(*allocLen = k) );
	else if ( *allocLen < k )
		*a = (int*)realloc(*a , sizeof(int)*(*allocLen = k));
}

static void box_partition_maps( 
  const int np ,
  const int my_p ,
  const Box* pbox,
  const int ghost ,
  int ** map_local_id ,
  int ** map_recv_pc ,
  int ** map_send_pc ,
  int ** map_send_id )
{
  const Box& my_box = pbox[my_p] ;

  const int my_ix = my_box[0][0] ;
  const int my_iy = my_box[1][0] ;
  const int my_iz = my_box[2][0] ;
  const int my_nx = my_box[0][1] - my_box[0][0] ;
  const int my_ny = my_box[1][1] - my_box[1][0] ;
  const int my_nz = my_box[2][1] - my_box[2][0] ;

  const int my_use_nx = 2 * ghost + my_nx ;
  const int my_use_ny = 2 * ghost + my_ny ;
  const int my_use_nz = 2 * ghost + my_nz ;

  const int id_length = my_use_nx * my_use_ny * my_use_nz ;

  int * local_id  = (int *) malloc( id_length * sizeof(int) );
  int * recv_pc   = (int *) malloc( ( np + 1 ) * sizeof(int) );
  int * send_pc   = (int *) malloc( ( np + 1 ) * sizeof(int) );

  int * send_id  = NULL ;
  int   send_id_size = 0 ;

  int iLocal , iSend ;
  int i ;

  Box my_use_box;

  my_use_box[0][0] = my_box[0][0] - ghost ;
  my_use_box[0][1] = my_box[0][1] + ghost ;
  my_use_box[1][0] = my_box[1][0] - ghost ;
  my_use_box[1][1] = my_box[1][1] + ghost ;
  my_use_box[2][0] = my_box[2][0] - ghost ;
  my_use_box[2][1] = my_box[2][1] + ghost ;

  for ( i = 0 ; i < id_length ; ++i ) { local_id[i] = -1 ; }

  iSend = 0 ;
  iLocal = 0 ;

  /* The vector space is partitioned by processors */

  for ( i = 0 ; i < np ; ++i ) {
    const int ip = ( i + my_p ) % np ;
    recv_pc[i] = iLocal ;
    send_pc[i] = iSend ;

    if ( ! box_disjoint( my_use_box , pbox[ip] ) ) {
      const int p_ix = pbox[ip][0][0] ;
      const int p_iy = pbox[ip][1][0] ;
      const int p_iz = pbox[ip][2][0] ;
      const int p_ex = pbox[ip][0][1] ;
      const int p_ey = pbox[ip][1][1] ;
      const int p_ez = pbox[ip][2][1] ;

      int local_x , local_y , local_z ;

      /* Run the span of global cells that my processor uses */

      for ( local_z = -ghost ; local_z < my_nz + ghost ; ++local_z ) {
      for ( local_y = -ghost ; local_y < my_ny + ghost ; ++local_y ) {
      for ( local_x = -ghost ; local_x < my_nx + ghost ; ++local_x ) {

        const int global_z = local_z + my_iz ;
        const int global_y = local_y + my_iy ;
        const int global_x = local_x + my_ix ;

        const int entry = 
          box_map_local_entry(my_box,ghost,local_x,local_y,local_z);

        if ( entry < 0 ) { abort(); }

        if ( p_iz <= global_z && global_z < p_ez &&
             p_iy <= global_y && global_y < p_ey &&
             p_ix <= global_x && global_x < p_ex ) {

          /* This ordinal is owned by processor 'ip' */

          local_id[ entry ] = iLocal++ ;

#if defined(DEBUG_PRINT)
if ( my_p != ip ) {
  fprintf(stdout,"  (%d,%d,%d) : P%d recv at local %d from P%d\n",
                  global_x,global_y,global_z,my_p,local_id[entry],ip);
  fflush(stdout);
}
#endif
        }

        /* If in my ownership and used by the other processor */
        if ( my_p != ip &&
             /* In my ownership: */
             ( 0 <= local_z && local_z < my_nz &&
               0 <= local_y && local_y < my_ny &&
               0 <= local_x && local_x < my_nx ) &&
             /* In other processors usage: */
             ( p_iz - ghost <= global_z && global_z < p_ez + ghost &&
               p_iy - ghost <= global_y && global_y < p_ey + ghost &&
               p_ix - ghost <= global_x && global_x < p_ex + ghost ) ) {

          resize_int( & send_id , & send_id_size , (iSend + 1) );
          send_id[ iSend ] = local_id[ entry ] ;
          ++iSend ;

#if defined(DEBUG_PRINT)
{
  fprintf(stdout,"  (%d,%d,%d) : P%d send at local %d to P%d\n",
                  global_x,global_y,global_z,my_p,local_id[entry],ip);
  fflush(stdout);
}
#endif
        }
      }
    }
    }
    }
  }
  recv_pc[np] = iLocal ;
  send_pc[np] = iSend ;

  *map_local_id  = local_id ;
  *map_recv_pc   = recv_pc ;
  *map_send_pc   = send_pc ;
  *map_send_id   = send_id ;
}

void box_partition_rcb( const int np , 
                        const int my_p ,
                        const Box& root_box,
                        const int ghost ,
                        Box** pbox,
                        int ** map_local_id ,
                        int ** map_recv_pc ,
                        int ** map_send_pc ,
                        int ** map_send_id )
{
  *pbox = new Box[ np ];

  box_partition( 0 , np , 2 , root_box , *pbox );

  box_partition_maps( np , my_p , *pbox , ghost ,
                      map_local_id , map_recv_pc , 
                      map_send_pc , map_send_id );
}

/*--------------------------------------------------------------------*/

#ifdef UNIT_TEST

static int box_contain( const Box& a , const Box& b )
{
  return a[0][0] <= b[0][0] && b[0][1] <= a[0][1] &&
         a[1][0] <= b[1][0] && b[1][1] <= a[1][1] &&
         a[2][0] <= b[2][0] && b[2][1] <= a[2][1] ;
}

static void box_print( FILE * fp , const Box& a )
{
  fprintf(fp,"{ [ %d , %d ) , [ %d , %d ) , [ %d , %d ) }",
                a[0][0] , a[0][1] ,  
                a[1][0] , a[1][1] ,  
                a[2][0] , a[2][1] );
}

static void test_box( const Box& box , const int np )
{
  const int ncell_box = box[0][1] * box[1][1] * box[2][1] ;
  int ncell_total = 0 ;
  int ncell_min = ncell_box ;
  int ncell_max = 0 ;
  std::vector<Box> pbox(np);
  int i , j ;

  box_partition( 0 , np , 2 , box , &pbox[0] );

  for ( i = 0 ; i < np ; ++i ) {
    const int ncell = ( pbox[i][0][1] - pbox[i][0][0] ) *
                      ( pbox[i][1][1] - pbox[i][1][0] ) *
                      ( pbox[i][2][1] - pbox[i][2][0] );

    if ( ! box_contain( box , pbox[i] ) ) {
      fprintf(stdout,"  OUT OF BOUNDS pbox[%d/%d] = ",i,np);
      box_print(stdout,pbox[i]);
      fprintf(stdout,"\n");
      abort();
    }

    for ( j = i + 1 ; j < np ; ++j ) {
      if ( ! box_disjoint( pbox[i] , pbox[j] ) ) {
        fprintf(stdout,"  NOT DISJOINT pbox[%d/%d] = ",i,np);
        box_print(stdout, pbox[i]);
        fprintf(stdout,"\n");
        fprintf(stdout,"               pbox[%d/%d] = ",j,np);
        box_print(stdout, pbox[j]);
        fprintf(stdout,"\n");
        abort();
      }
    }
    ncell_total += ncell ;

    if ( ncell_max < ncell ) { ncell_max = ncell ; }
    if ( ncell < ncell_min ) { ncell_min = ncell ; }
  }

  if ( ncell_total != ncell_box ) {
    fprintf(stdout,"  WRONG CELL COUNT NP = %d\n",np);
    abort();
  }
  fprintf(stdout,"NP = %d, total = %d, avg = %d, min = %d, max = %d\n",
          np,ncell_box,ncell_box/np,ncell_min,ncell_max);
}

/*--------------------------------------------------------------------*/

#endif


