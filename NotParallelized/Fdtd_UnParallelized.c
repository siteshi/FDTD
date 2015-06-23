/* FD1D_1.1.c.  1D FDTD simulation in free space */

# include <math.h>
# include <stdlib.h>
# include <stdio.h>

#define KE  200         /* KE is the number of cells to be used */

int main ()
{
  float ex[KE], hy[KE];
  int n, k, kc, ke, NSTEPS;
  float T;
  float t0, spread, pulse;
  FILE *fpE,*fpH, *fopen();

  /* Initialize */
  for ( k=0; k < KE; k++ ) {
    ex[k] = 0;
    hy[k] = 0.;
  }

  kc = KE/2;    /* Center of the problem space */
 

  t0 = 40.0;    /* Center of the incident pulse */
  spread = 12;  /* Width of the incident pulse */
  T = 0;
  NSTEPS = 1;



  while ( NSTEPS > 0 ) {
    printf( "NSTEPS --> ");     /* NSTEPS is the number of times the */
    scanf("%d", &NSTEPS);       /* main loop has executed */
    printf("%d \n", NSTEPS);
    n = 0;
  fpE = fopen( "Ex", "w");
  fpH = fopen( "Hy", "w");
    for (n = 1; n <= NSTEPS; n++) {
      T = T + 1;        /* T keeps track of the total number */
                        /* of times the main loop is executed */
      /* Main FDTD Loop */

        /* Calculate the Ex field */
        ex[0] = ex[0] + 0.5*(0 - hy[0]);
        for ( k = 1; k < KE; k++ ) {
          ex[k] = ex[k] + 0.5*(hy[k-1] - hy[k]);
        }

        /* Put a Gaussian pulse in the middle */
        pulse = exp(-0.5*(pow( (t0-T)/spread, 2.0 )));
         kc = 125;
        ex[kc] = pulse;
 	//ex[kc + 20] = pulse;
	
	
        printf("%5f  %6f\n", t0-T, ex[kc]);

        /* Calculate the Hy field */
        for ( k = 0; k < KE-1; k++ ) {
          hy[k] = hy[k] + .5*(ex[k] - ex[k+1] );
        }
	 hy[KE-1] = hy[KE-1] + .5*(ex[KE-1] - 0 );

	/* Write the E field out to a file "Ex" */
	//fprintf( fpE, "  Start of Time %6d \n", n);
	for ( k=0; k<=KE; k++ )
      { fprintf( fpE, "%f,", ex[k]); }
	fprintf( fpE, "\n");


	 /* Write the H field out to a file "Hy" */
	//fprintf( fpH, "  Start of Time %6d \n", n);
	for ( k=0; k<=KE; k++ )
      { // fprintf( fpH, "  %6.2f \n", hy[k]); 
	 fprintf( fpH, "%f,", hy[k]);
	}
	fprintf( fpH, "\n");

    }

    /* End of the Main FDTD Loop */

      /* At the end of the calculation, print out
            the Ex and Hy fields */
      for ( k=0; k<= KE; k++ )
      { printf( "%3d    %6f   %6f\n", k, ex[k], hy[k]); }

      
    
      fclose(fpE);
      fclose(fpH);

      printf( "T = %5.0f\n", T);

  }
}

