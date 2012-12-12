/*
 * Gregory Petropoulos
 *
 * This is my Dilute Ising Model program for the final
 *
 * To compile:  mpicc -g -Wall -std=c99 -o DIsing dising.c
 * To run:  mpiexec -n 1 ./DIsing <input file name> <partition> <beta> <mu> <run iterations> <measurment iterations> <save iterations>
 *          <> -> mandatory
 *          [] -> optional
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "mpi.h"
#include "lib/ran3.h"

// -----------------------------------------------------------------
// global variables --> perhaps make some of these local?
int nprocs;
int my_rank;
int field_width;                                        /* local board size                       */
int field_height; 
int local_width;                                        /* local size of data                     */
int local_height;
int width;                                              /* The total dimension of the field       */
int height;
int ncols;
int nrows;
int *field;                                             /* The local data fields                  */
char header[100];
MPI_Datatype row;                                     /* makes row data type                    */
MPI_Datatype col;                                     /* makes col data type                    */
MPI_Datatype etype, filetype, contig;                 /* derrived data types for IO             */
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Ends the program on an error and prints message to node 0
void cleanup (const char *message) {
  if (my_rank == 0) printf("%s\n",message);
  MPI_Type_free(&row);
  MPI_Type_free(&col);
  MPI_Finalize();                                       /* kills mpi                              */
  exit(0);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// reads in file
void fileread (char *filename, char *partition, int *offset) {
  if (strcmp(partition, "slice") == 0) {                /* determines the data decomposition      */
    ncols = 1;
    nrows = nprocs;
  }
  else if (strcmp(partition, "checkerboard") == 0) {
    ncols = sqrt(nprocs);
    nrows = sqrt(nprocs);
  }
  else {
    ncols = 1;
    nrows = 1;
  }

  MPI_File fh;                                          /* sets up MPI input                      */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  int newline_count = 0;
  char read;
  *offset = 0;

  do {                                                  /* reads in the header                    */
    MPI_File_read_at_all(fh, *offset, &read, 1, MPI_CHAR, MPI_STATUS_IGNORE);
    header[*offset]=read;
    if (read == '\n') newline_count++;
    (*offset)++;
    if (*offset == 100) cleanup("Error:  Header exceeds 100 characters, check file or recompile code");
  } while (newline_count < 3);
    
  char head[10];                                        /* parses the header and error            */
  int depth;                                            /*   checks                               */
  int rv               =  sscanf(header, "%6s\n%i %i\n%i\n", head, &width, &height, &depth);
  if (rv != 4)            cleanup("Error: The file did not have a valid PGM header");
  if (my_rank==0)         printf( "%s: %s %i %i %i\n", filename, head, width, height, depth );
  if (strcmp(head, "P5")) cleanup("Error: PGM file is not a valid P5 pixmap");
  if( depth != 255 )      cleanup("Error: PGM file requires depth=255");
  if (width % ncols)      cleanup("Error: pixel width cannot be divided evenly into cols");
  if (height % nrows)     cleanup("Error: %i pixel height cannot be divided into %i rows\n");

  local_width = width / ncols;                          /* determines the size of                 */
  local_height = height / nrows;                        /* the local data                         */
  field_width = local_width + 2;                        /* creates an array with room for         */
  field_height = local_height + 2;                      /* ghosts and boarders                    */
  // playing fields
  field = (int *)malloc(field_width * field_height * sizeof(int));
  // array for the file read
  char *temp=(char *)malloc( local_width * local_height * sizeof(char));
  MPI_Aint extent;                                      /* declares the extent                    */
  MPI_Offset disp = *offset;                            /* the initial displacement of            */
                                                        /*   the header                           */
  // this needs to be added to the displacement so each processor starts reading from
  // the right point of the file
  disp += (my_rank / ncols) * width * local_height + (my_rank % ncols) * local_width;
  
  etype = MPI_CHAR;                                     /* sets the etype to MPI_CHAR             */
  MPI_Type_contiguous(local_width, etype, &contig);     /*                                        */
  extent = width * sizeof(char);                        /* total size of repeatable unit          */
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype);                           /* makes the filetype derrived data       */
  // reads in the file
  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, temp, local_width*local_height, MPI_CHAR, MPI_STATUS_IGNORE);

  int x, y, b, ll; 
  for (y = 0; y < field_height; y++ ) {                 /* loops through the field                */
    for(x = 0; x < field_width; x++ ) {
      ll = (y * field_width + x);                       /* puts zeros at the boarders             */
      if ((x == 0) || (y == 0) || (x == field_width - 1) || (y == field_height - 1)) {
        field[ll] = 0;
      }
      else {
        b = (int)temp[x-1+(y-1)*local_width];           /* finds data from read                   */
        field[ll] = b - 2;
      }
    }
  }
  MPI_Type_free(&filetype);
  MPI_File_close(&fh);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// writes the game board out to file
void filewrite (char *in_file, int iteration, int offset, int total_iterations) {
  char filename[1000];
  sprintf(filename, "Configs/%0*d_", (int)log10((double) total_iterations)+1, iteration); 
  strcat(filename, in_file);

  MPI_File fh;                                          /* sets up MPI input                      */
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  if (my_rank == 0) {                                   /* writes the header                      */
    MPI_File_write(fh,header,offset,MPI_CHAR,MPI_STATUS_IGNORE); 
  }

  char *temp=(char *)malloc( local_width * local_height * sizeof(char));

  int x, y, ll;
  for (y = 0; y < field_height; y++ ) {                 /* loops through the field                */
    for(x = 0; x < field_width; x++ ) {
      if ((x != 0) && (y != 0) && (x != field_width - 1) && (y != field_height - 1)) {
        ll = (y * field_width + x);                     /* puts zeros at the boarders             */
        temp[x-1+(y-1)*local_width] = (char)(field[ll] + 2);
        /*b = field_pointer[ll];*/
        /*b = (b==0)?0xFF:0;                              [> black = bugs; other = no bug           <]*/
        /*temp[x-1+(y-1)*local_width] = (char)b;*/
      }
    }
  }

  MPI_Aint extent;                                      /* declares the extent                    */
  MPI_Offset disp = offset;                             /* the initial displacement of            */

  // this needs to be added to the displacement so each processor starts reading from
  // the right point of the file
  disp += (my_rank / ncols) * width * local_height + (my_rank % ncols) * local_width;

  etype = MPI_CHAR;                                     /* sets the etype to MPI_CHAR             */
  MPI_Type_contiguous(local_width, etype, &contig);     /*                                        */
  extent = width * sizeof(char);                        /* total size of repeatable unit          */
  MPI_Type_create_resized(contig, 0, extent, &filetype); 
  MPI_Type_commit(&filetype);                           /* makes the filetype derrived data       */
  // writes the file
  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, temp, local_width*local_height, MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_Type_free(&filetype);
  MPI_File_close(&fh);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// computes <m>, <e>
void measure (int iteration, double mu) {
  int local_m = 0;
  int global_m = 0;
  int local_m2 = 0;
  int global_m2 = 0;
  int local_e = 0;
  int global_e = 0;
  int i, j, top, right;

  for (j = 0; j < local_height; j++) {
    for (i = 0; i < local_width; i++) {
      local_m += field[(j+1)*(field_width)+(i+1)];
      local_m2 += field[(j+1)*(field_width)+(i+1)]*field[(j+1)*(field_width)+(i+1)];
      top = field[(j)*(field_width)+(i+1)];
      right = field[(j+1)*(field_width)+(i+2)];
      local_e += (field[(j+1)*(field_width)+(i+1)]*(right+top));
    }
  }
  
  if (nprocs > 1) {
    MPI_Reduce(&local_m, &global_m, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_e, &global_e, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_m2, &global_m2, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  else {
    global_m = local_m;
    global_e = local_e;
    global_m2 = local_m2;
  }
  global_e += mu*local_m2;

  if (my_rank == 0) {printf("MEASUREMENT for sweep %i: %e %e\n", iteration, (double)global_m/(double)(width*height) , (double)global_e/(double)(width*height));}
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// swaps ghost rows
void summonspectre(MPI_Request* send, MPI_Request* recv) {

  /* POINTER LOCATIONS *
   *                   *
   *  |a|b| ... |c|d|  *
   *  |e|              *
   *  |.|              *
   *  |.|              *
   *  |.|              *
   *  |f|              *
   *  |g|              *
   *                   *
   *********************/

  int *p_a = field;                                   /* pointers for communication             */
  int *p_b = field + 1;
  int *p_c = field + local_width;
  int *p_d = field + local_width + 1;
  int *p_e = field + field_width;
  int *p_f = field + field_width * local_height;
  int *p_g = field + field_width * (local_height + 1);

  int left = (my_rank%ncols==0?ncols:0)+my_rank-1;
  int right = (my_rank%ncols==ncols-1?-ncols:0)+my_rank+1;
  int up = (my_rank/ncols==0?ncols*nrows:0)+my_rank-ncols;
  int down = (my_rank/ncols==nrows-1?-ncols*nrows:0)+my_rank+ncols;

  /*int x = my_rank%ncols;*/
  /*int y = my_rank/ncols);*/

  /*int up = (x + ((y+1)%nrows)*ncols);*/
  /*int right = ((x+1)%ncols+y*ncols);*/
  /*int down = (x + ((y+nrows-1)%nrows)*ncols);*/
  /*int left = ((x+ncols-1)%ncols+y*ncols);*/

  MPI_Isend(p_b, 1, col, left, 0, MPI_COMM_WORLD, &send[0]);
  MPI_Irecv(p_d, 1, col, right, 0, MPI_COMM_WORLD, &recv[0]);
  MPI_Isend(p_c, 1, col, right, 0, MPI_COMM_WORLD, &send[1]);
  MPI_Irecv(p_a, 1, col, left, 0, MPI_COMM_WORLD, &recv[1]);

  MPI_Waitall(4, recv, MPI_STATUS_IGNORE);

  MPI_Isend(p_e, 1, row, up, 0, MPI_COMM_WORLD, &send[2]);
  MPI_Irecv(p_g, 1, row, down, 0, MPI_COMM_WORLD, &recv[2]);
  MPI_Isend(p_f, 1, row, down, 0, MPI_COMM_WORLD, &send[3]);
  MPI_Irecv(p_a, 1, row, up, 0, MPI_COMM_WORLD, &recv[3]);
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// generates the probabilities needed for the lookup tables
void lookup (double beta, double mu, double *p_up, double *p_dn, double *p_no) {
  int i;
  double e_up, e_dn, e_no, norm, sum;

  for (i = 0; i<9; i++) {
    e_up = exp(beta*((i-4)+mu));
    e_dn = exp(beta*(-1*(i-4)+mu));
    e_no = 1;
    norm = e_up + e_dn + e_no;

    p_up[i] = e_up / norm;
    p_dn[i] = e_dn / norm;
    p_no[i] = e_no / norm;

    sum = p_up[i] + p_dn[i] + p_no[i];

    if (my_rank == 0) {printf("LOOKUP for nn sum %i:  %e + %e + %e = %e\n", i-4, p_up[i], p_dn[i], p_no[i], sum);}
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// updates the spins
void update (const double *p_up, const double *p_dn, const double *p_no) {
  int i, x, y, xb, yb, neighbor;
  int seed = 1;
  double ran; 
  for (i = 0; i < local_height*local_width; i++) {    /* loops through the local data           */
    ran = rand() / (double)RAND_MAX;
    /*ran = ran3(&seed);*/
    x = (int)(ran*(double)local_width);
    ran = rand() / (double)RAND_MAX;
    /*ran = ran3(&seed);*/
    y = (int)(ran*(double)local_height);

    yb = y + 1;                                       /* shifts needed becuase the board        */
    xb = x + 1;                                       /*   is bigger due to ghost rows          */

    neighbor  = field[(yb-1)*field_width+xb];         /* top                                    */
    neighbor += field[(yb+1)*field_width+xb];         /* bottom                                 */
    neighbor += field[(yb)*field_width+xb-1];         /* left                                   */
    neighbor += field[(yb)*field_width+xb+1];         /* right                                  */

    ran = rand() / (double)RAND_MAX;
    /*ran = ran3(&seed);*/

    if ((0 <= ran) && (ran < p_up[neighbor+4])) {field[yb*field_width+xb] = 1;}
    else if ((p_up[neighbor+4] <= ran) && (ran < p_dn[neighbor+4]+p_up[neighbor+4])) {field[yb*field_width+xb] = -1;}
    else {field[yb*field_width+xb] = 0;}
  }
}
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// the main program
int main (int argc, char *argv[]) {
  char in_file[1000];
  char partition[100];
  int r_iterations, m_interval, w_interval, offset, i;
  double beta, mu;
  double p_up[9], p_dn[9], p_no[9];
  MPI_Request send[4], recv[4];
  int seed = -100-my_rank;
  ran3(&seed);

  MPI_Init(&argc, &argv);                               /* start up MPI                           */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);               /* get the number of processes            */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);              /* get my rank among all the processes    */


  for (i = 0; i < 4; i++) {                             /* initializes the send and recv Request  */
    send[i] = MPI_REQUEST_NULL;                         /*  variables to NULL                     */
    recv[i] = MPI_REQUEST_NULL;
  }
  if (argc < 8) {                                       /* parses command line arguments          */
    cleanup("Error:  Too few arguments");
  }
  else if (argc == 8) {
    strcpy(in_file,   argv[1]);
    strcpy(partition, argv[2]);
    if (!((strcmp(partition, "slice") == 0) || (strcmp(partition, "checkerboard") == 0) || (strcmp(partition, "none") == 0))) {
      cleanup("Error:  Incorrect partition option.  Enter 'slice', 'checkerboard', or 'none'");
    }
    beta = strtod(argv[3], NULL);
    mu = strtod(argv[4], NULL);
    r_iterations = atoi(argv[5]);
    m_interval = atoi(argv[6]);
    w_interval = atoi(argv[7]);
  }
  else if (argc > 8) {
    cleanup("Error:  Too many arguments");
  }

  fileread(in_file, partition, &offset);                /* function call to read the file in      */
  lookup(beta, mu, p_up, p_dn, p_no);

  MPI_Type_vector(field_height, 1, field_width, MPI_INT, &col);
  MPI_Type_commit(&col);

  MPI_Type_vector(field_width, 1, 1, MPI_INT, &row);
  MPI_Type_commit(&row);

  for (i = 0; i < r_iterations; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    summonspectre(send, recv);                          /* exchanges ghost fields                 */
    if (m_interval > 0) {
      if (i%m_interval == 0) {
        measure(i, mu);
      }
    }

    if (w_interval > 0) {
      if (i%w_interval == 0) {
        filewrite(in_file, i, offset, r_iterations);
      }
    }

    MPI_Waitall(4, recv, MPI_STATUS_IGNORE);            /* waits on sends                         */
    update(p_up, p_dn, p_no);                           /* updates the spins                      */
    MPI_Waitall(4, send, MPI_STATUS_IGNORE);            /* waits on receives                      */
  }
  cleanup("Run complete!");                    /* closes the program                     */
  return 0;
}
