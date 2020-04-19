#include "sph.h"

#define CAP_INCREMENT 4096

bool verbose = false, restart = false;
int size, rank, west, east;
real_t subdomain[2];
int_t n_field = 0,
        n_virt = 0,
        n_mirror = 0,
        n_pairs = 0,
        n_global_field = 0;

/* Timing */
double
    t_start, t_end,
    t_generate = 0.0,
    t_border = 0.0,
    t_timestep = 0.0,
    t_migrate = 0.0,
    t_io = 0.0,
    t_find_neighbors = 0.0;

/* Run parameters */
int_t
    min_iteration = MIN_ITERATION_DEFAULT,
    max_iteration = MAX_ITERATION_DEFAULT,
    checkpoint_frequency = CHECKPOINT_FREQUENCY_DEFAULT;

particle_t *list;                   // Flat list of particle structs
int_t n_capacity = CAP_INCREMENT;   // Initial list capacity, grows

pair_t *pairs;
int_t n_pair_cap = CAP_INCREMENT;

bucket_t** buckets;


void print_timing(char* full_string, char* short_string, double value) {
    if (rank == 0) {
        printf(full_string, value);
        for (int r = 1; r < size; ++r) {
            size_t remote_string_size = strlen(short_string);
            char remote_string[remote_string_size];
            double remote_value = 0.0;
            MPI_Recv(remote_string, remote_string_size+1, MPI_CHAR, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&remote_value, 1, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf(remote_string, remote_value);
        }
        printf("\n");
    }
    else {
        size_t string_size = strlen(short_string);
        MPI_Send(short_string, string_size+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}

void
ext_force ( void )
{
    #pragma omp parallel for
    for ( int_t k=0; k<n_field; k++ )
        EXDVXDT(k,1) = -9.81;
}


void
int_force ( void )
{

    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int_t k = 0; k < (n_field + n_virt + n_mirror); k++)
            INDVXDT(k, 0) = INDVXDT(k, 1) = 0.0;

        #pragma omp for
        for (int_t k = 0; k < n_field; k++)
            if (INTER(k) < free_surface)
                RHO(k) = density;

        #pragma omp for
        for (int_t k = 0; k < (n_field + n_virt + n_mirror); k++)
            P(k) = sos * sos * density * ((pow(RHO(k) / density, 7.0) - 1.0) / 7.0);

        // All the pairwise interactions
        #pragma omp for
        for (int x = 0; x < N_BUCKETS_X; ++x) {
            for (int y = 0; y < N_BUCKETS_Y; ++y) {
                pair_t *p = &pairs[BID(x, y)];
                while (p != NULL && p->ip != NULL) {
                    int_t
                            i = p->i,
                            j = p->j;
                    real_t hx, hy;

                    // i acts on j
                    hx = -(P(i) / pow(RHO(i), 2) + P(j) / pow(RHO(j), 2)) * p->dwdx[0];
                    hy = -(P(i) / pow(RHO(i), 2) + P(j) / pow(RHO(j), 2)) * p->dwdx[1];

                    INDVXDT(i, 0) += M(j) * hx;

                    INDVXDT(i, 1) += M(j) * hy;

                    if (p->jp->bucket_x == x && p->jp->bucket_y == y) {
                        // j acts on i, reverse sign because dwdx is X(i)-X(j)
                        hx = -(P(j) / pow(RHO(j), 2) + P(i) / pow(RHO(i), 2)) * (-p->dwdx[0]);
                        hy = -(P(j) / pow(RHO(j), 2) + P(i) / pow(RHO(i), 2)) * (-p->dwdx[1]);

                        INDVXDT(j, 0) += M(i) * hx;

                        INDVXDT(j, 1) += M(i) * hy;
                    }


                    p = p->next;

                }
            }
        }
    }

}


void
correction ( void ) {

    #pragma omp parallel
    {
        #pragma omp for
        for (int x = 0; x < N_BUCKETS_X; ++x) {
            for (int y = 0; y < N_BUCKETS_Y; ++y) {
                pair_t *p = &pairs[BID(x, y)];
                while (p != NULL && p->ip != NULL) {
                    int_t
                            i = p->i,
                            j = p->j;
                    real_t drho;
                    drho = RHO(i) - RHO(j);
                    AVRHO(i) -= drho * p->w / WSUM(i);
                    if (p->jp->bucket_x == x && p->jp->bucket_y == y) {
                        drho = RHO(j) - RHO(i);
                        AVRHO(j) -= drho * p->w / WSUM(j);
                    }
                    p = p->next;
                }
            }
        }

        #pragma omp for
        for (int_t k = 0; k < (n_field + n_virt + n_mirror); k++) {
            if (TYPE(k) < 0 && INTER(k) < 10)
                RHO(k) = density;
            else
                RHO(k) += 0.5 * AVRHO(k);
        }
    }
}


void
cont_density ( void )
{
    for ( int_t k=0; k<(n_field+n_virt+n_mirror); k++ )
        DRHODT(k) = 0.0;


    #pragma omp parallel for
    for (int x = 0; x < N_BUCKETS_X; ++x) {
        for (int y = 0; y < N_BUCKETS_Y; ++y) {
            pair_t* p = &pairs[BID(x, y)];
            while (p != NULL && p->ip != NULL) {
                int_t
                    i = p->i,
                    j = p->j;
                real_t vcc;

                vcc = (VX(i)-VX(j))*p->dwdx[0] +
                    (VY(i)-VY(j))*p->dwdx[1];

                DRHODT(i) += RHO(i) * (M(j)/RHO(j)) * vcc;

                // Reverse sign of dwdx because it is calc. with X(i)-X(j)
                if (p->jp->bucket_x == x && p->jp->bucket_y == y) {
                    vcc = (VX(j)-VX(i))*(-p->dwdx[0]) +
                        (VY(j)-VY(i))*(-p->dwdx[1]);

                    DRHODT(j) += RHO(j) * (M(i)/RHO(i)) * vcc;
                }


                p = p->next;
            }

        }
    }

    #pragma omp parallel for
    for ( int_t k=0; k<(n_field+n_virt+n_mirror); k++ )
    {
        RHO(k) += 0.5 * dt * DRHODT(k);
    }
}


void
kernel ( void )
{
    static real_t factor = 7.0 / (478.0 * M_PI * H * H);

    #pragma omp parallel for
    for (int x = 0; x < N_BUCKETS_X; ++x) {
        for (int y = 0; y < N_BUCKETS_Y; ++y) {
            pair_t* p = &pairs[BID(x, y)];
            while(p != NULL && p->ip != NULL) {
                real_t q = p->q;
                real_t dx[2] = {p->ip->x[0] - p->jp->x[0], p->ip->x[1] - p->jp->x[1]};

                if ( q == 0.0 )
                {
                    p->w = factor * (
                                     pow((3-q),5) - 6*pow((2-q),5) + 15*pow((1-q),5)
                                     );
                    p->dwdx[0] = p->dwdx[1] = 0.0;
                }
                else if ( q>0.0 && q<=1.0 )
                {
                    p->w = factor * (
                                     pow((3-q),5) - 6*pow((2-q),5) + 15*pow((1-q),5)
                                     );
                    p->dwdx[0] = (factor/pow(H,2)) *
                        (-120+120*pow(q,2)-50*pow(q,3))*dx[0];
                    p->dwdx[1] = (factor/pow(H,2)) *
                        (-120+120*pow(q,2)-50*pow(q,3))*dx[1];
                }
                else if ( q>1.0 && q<=2.0 )
                {
                    p->w = factor * ( pow(3-q,5) - 6*pow(2-q,5));
                    p->dwdx[0] = (factor/H) *
                        ((-5)*pow((3-q),4)+30*pow((2-q),4))*(dx[0]/p->r);
                    p->dwdx[1] = (factor/H) *
                        ((-5)*pow((3-q),4)+30*pow((2-q),4))*(dx[1]/p->r);
                }
                else if ( q>2.0 && q<=3.0 )
                {
                    p->w = factor * pow(3-q,5);
                    p->dwdx[0] = (factor/H) * ((-5)*pow((3-q),4))*(dx[0]/p->r);
                    p->dwdx[1] = (factor/H) * ((-5)*pow((3-q),4))*(dx[1]/p->r);
                }
                else
                {
                    p->w = 0.0;
                    p->dwdx[0] = p->dwdx[1] = 0.0;
                }
                WSUM(p->i) += p->w;
                if (p->jp->bucket_x == x && p->jp->bucket_y == y) {
                    WSUM(p->j) += p->w;
                }

                p = p->next;
            }
        }
    }
}

void find_neighbors_buckets_ws( void ) {
    int_t n_total = n_field + n_virt + n_mirror;
    n_pairs = 0;


    #pragma omp parallel for
    for ( int_t k=0; k<n_total; k++ )
        INTER(k) = WSUM(k) = AVRHO(k) = 0;


    if (buckets == NULL) {
        fprintf(stderr, "4Not enough memory! Line %d\n", __LINE__);
        exit(1);
    }

#ifdef FILL_BUCKETS_LOCK
    omp_lock_t lock[N_BUCKETS_X*N_BUCKETS_Y];

    for (int i=0; i<N_BUCKETS_X*N_BUCKETS_Y; i++)
        omp_init_lock(&(lock[i]));
#endif //FILL_BUCKETS_LOCK

    #pragma omp parallel
    {
        /* Init buckets */
        #pragma omp for nowait
        for (int x = 0; x < N_BUCKETS_X; ++x) {
            for (int y = 0; y < N_BUCKETS_Y; ++y) {
                buckets[BID(x, y)]->particle = NULL;
                buckets[BID(x, y)]->next = NULL;
            }
        }

        /* Compute bucket_x and bucket_y for all particles */
        #pragma omp for
        for (int i = 0; i < n_total; ++i) {
            particle_t *particle = &list[i];

            int actual_x = MIN((int) (((particle->x[0] - subdomain[0])+RADIUS) / BUCKET_RADIUS) , N_BUCKETS_X-1);
            int actual_y = MIN((int) ((particle->x[1]+1.55*H) / BUCKET_RADIUS),N_BUCKETS_Y-1);

            particle->local_idx = i;
            particle->bucket_x = actual_x;
            particle->bucket_y = actual_y;
        }

        /* Fill buckets */
#ifdef FILL_BUCKETS_LOCK
        fill_buckets2(buckets, n_total, lock);
#else
        fill_buckets1(buckets, n_total);
#endif //FILL_BUCKETS_LOCK

        /* Create neighbors 2 */
        #pragma omp for
        for (int x = 0; x < N_BUCKETS_X; ++x) {
            for (int y = 0; y < N_BUCKETS_Y; ++y) {
                pair_t* base_pair = &pairs[BID(x, y)];
                bucket_t* base_bucket = buckets[BID(x, y)];
                particle_t* base_particle;
                if (base_bucket != NULL) { //kan kanskje fjernes
                    base_particle = base_bucket->particle;
                    if (base_particle == NULL) {
                        continue;
                    }
                }
                else
                    continue;

                while(base_bucket != NULL && base_bucket->particle != NULL) {
                    base_particle = base_bucket->particle;

                    int noffsets[8][2] = { //kan kanskje flyttes opp
                                                 {-1, 1}, // Top left
                                                 {0, 1}, // Top center
                                                 {1, 1},  // Top right
                                                 {1, 0},  // Right
                                                 {1, -1}, // Bottom right
                                                 {0, -1}, // Bottom center
                                                 {-1, -1}, // Bottom left
                                                 {-1, 0},  // Left
                    };

                    for (int i = 0; i < 9; ++i) {
                        bucket_t* current_bucket;
                        if (i == 0) {
                            current_bucket = base_bucket->next;
                        } else {
                            int bx = x+noffsets[i-1][0];
                            int by = y+noffsets[i-1][1];
                            if (bx >= N_BUCKETS_X || bx < 0 ||
                                by >= N_BUCKETS_Y || by < 0) continue;
                            current_bucket =
                                buckets[BID(bx, by)];
                        }
                        while (current_bucket != NULL && current_bucket->particle != NULL) {
                            particle_t* current_particle = current_bucket->particle;
                            double distance_squared =
                                pow(base_particle->x[0] - current_particle->x[0], 2) +
                                pow(base_particle->x[1] - current_particle->x[1], 2);
                            if (distance_squared <= RADIUS*RADIUS) {

                                if (base_pair->ip == NULL) {
                                    base_pair->i = base_particle->local_idx;
                                    base_pair->j = current_particle->local_idx;
                                    base_pair->ip = base_particle;
                                    base_pair->jp = current_particle;
                                    base_pair->r = sqrt(distance_squared);
                                    base_pair->q = base_pair->r / H;
                                    base_pair->w = 0.0;
                                    base_pair->dwdx[0] = base_pair->dwdx[1] = 0.0;
                                } else {
                                    base_pair->next = (pair_t*)malloc(sizeof(pair_t));
                                    base_pair->next->i = base_particle->local_idx;
                                    base_pair->next->j = current_particle->local_idx;
                                    base_pair->next->ip = base_particle;
                                    base_pair->next->jp = current_particle;
                                    base_pair->next->r = sqrt(distance_squared);
                                    base_pair->next->q = base_pair->next->r / H;
                                    base_pair->next->w = 0.0;
                                    base_pair->next->dwdx[0] = base_pair->next->dwdx[1] = 0.0;
                                    base_pair->next->next = NULL;
                                    base_pair = base_pair->next;
                                }

                                base_particle->interactions++;
                                if (i == 0) {
                                    current_particle->interactions++;
                                }

                            }
                            current_bucket = current_bucket->next;
                        }
                    }
                    base_bucket = base_bucket->next;
                }
            }
        }

        /* Free buckets */
        #pragma omp for
        for (int x = 0; x < N_BUCKETS_X; ++x) {
            for (int y = 0; y < N_BUCKETS_Y; ++y) {
                buckets[BID(x, y)]->particle = NULL;
                if (buckets[BID(x, y)]->next != NULL) {
                    free_buckets(buckets[BID(x, y)]->next);
                }
            }
        }
    }
#ifdef FILL_BUCKETS_LOCK
    for (int i=0; i<N_BUCKETS_X*N_BUCKETS_Y; i++)
        omp_destroy_lock(&(lock[i]));
#endif //FILL_BUCKETS_LOCK

}

static inline void fill_buckets1(bucket_t** buckets, int n_total) {
    #pragma omp for
    for (int x = 0; x < N_BUCKETS_X; ++x) {
        for (int y = 0; y < N_BUCKETS_Y; ++y) {
            bucket_t* bucket = buckets[BID(x, y)];

            for (int i = 0; i < n_total; ++i) {
                particle_t* particle = &list[i];

                if (particle->bucket_x != x || particle->bucket_y != y)  {
                    continue;
                }

                if(bucket->particle == NULL) {
                    bucket->particle = particle;
                } else {
                    bucket_t* new_bucket = (bucket_t*)malloc(sizeof(bucket_t));
                    new_bucket->particle = particle;
                    new_bucket->next = bucket;

                    buckets[BID(x, y)] = new_bucket;
                }
            }
        }
    }
}

#ifdef FILL_BUCKETS_LOCK
static inline void fill_buckets2(bucket_t** buckets, int n_total, omp_lock_t lock[N_BUCKETS_X*N_BUCKETS_Y]) {
    #pragma omp for
    for (int i = 0; i < n_total; ++i) {
        /* Give particle a local index */
        list[i].local_idx = i;

        /* Compute bucket coordinates */
        int bucket_x = list[i].bucket_x;
        int bucket_y = list[i].bucket_y;

        /*Lock possible critical section*/
        omp_set_lock(&(lock[BID(bucket_x,bucket_y)]));

        /* Pointer to relevant bucket */
        bucket_t *bucket = buckets[BID(bucket_x, bucket_y)];

        if (bucket->particle == NULL) {
            bucket->particle = &list[i];
            bucket->next = NULL;
        } else {
            bucket_t* new_bucket = (bucket_t*)malloc(sizeof(bucket_t));
            new_bucket->particle = &list[i];
            new_bucket->next = bucket;

            buckets[BID(bucket_x, bucket_y)] = new_bucket;
        }

        /*Unlock*/
        omp_unset_lock(&(lock[BID(bucket_x,bucket_y)]));

    }
}
#endif //FILL_BUCKETS_LOCK


void free_buckets(bucket_t* bucket) {
    if (bucket->next != NULL) {
        free_buckets(bucket->next);
    }
    free(bucket);
}



void
time_step ( int_t timestep )
{
    if ( timestep > 0 )
    {
        #pragma omp parallel for
        for ( int_t k=0; k<n_field; k++ )
        {
            // Vx, Vy cloned to vx_min vy_min for some reason
            // Omitting this until I see the purpose
            VX(k) += 0.5 * dt * DVX(k,0);
            VY(k) += 0.5 * dt * DVX(k,1);
            X(k) += dt * VX(k);
            Y(k) += dt * VY(k);
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < N_BUCKETS_X*N_BUCKETS_Y; ++i) {
        pairs[i].ip = NULL;
        pairs[i].jp = NULL;
        pairs[i].next = NULL;
        pairs[i].i = 0;
        pairs[i].j = 0;
    }


    MPI_Barrier(MPI_COMM_WORLD);
    double fn_start = MPI_Wtime();
    find_neighbors_buckets_ws();
    double fn_end = MPI_Wtime();
    t_find_neighbors += fn_end - fn_start;

    kernel();
    cont_density();
    if ( timestep > 0 )
        correction();
    int_force();
    ext_force();

    #pragma omp parallel for
    for ( int_t k=0; k<n_field; k++ )
    {
        DVX(k,0) = INDVXDT(k,0) + EXDVXDT(k,0); //+ ardvxdt
        DVX(k,1) = INDVXDT(k,1) + EXDVXDT(k,1); //+ ardvxdt
    }

    if ( timestep == 0 )
    {
        #pragma omp parallel for
        for ( int_t k=0; k<n_field; k++ )
        {
            // Calculate initial distributions
            RHO(k) += 0.5 * dt * DRHODT(k);
            VX(k) += 0.5 * dt * DVX(k,0);
            VY(k) += 0.5 * dt * DVX(k,1);
            X(k) += dt * VX(k);
            Y(k) += dt * VY(k);
        }
    }
    else
    {
        #pragma omp parallel for
        for ( int_t k=0; k<(n_field); k++ )
        {
            RHO(k) += 0.5 * dt * DRHODT(k);

            VX(k) += 0.5 * dt * DVX(k,0);
            VY(k) += 0.5 * dt * DVX(k,1);

            // Reflect velocity at boundaries
            if (Y(k) < 0.0 && VY(k) < 0.0)
                VY(k) = -VY(k);
            if (X(k) > B && VX(k) > 0.0)
                VX(k) = -VX(k);
            if (X(k) < 0.0 && VX(k) < 0.0)
                VX(k) = -VX(k);

            X(k) += dt * VX(k);
            Y(k) += dt * VY(k);

        }
    }

    #pragma omp parallel for
    for (int x = 0; x < N_BUCKETS_X; ++x) {
        for (int y = 0; y < N_BUCKETS_Y; ++y) {
            pair_t* p = &pairs[BID(x, y)];
            p->i = 0;
            p->j = 0;
            p->ip = NULL;
            p->jp = NULL;
            p = p->next;
            pairs[BID(x,y)].next = NULL;
            while(p != NULL) {
                pair_t* old = p;
                p = p->next;
                free(old);
            }
        }
    }
}


void
time_integration ( void )
{
    /* Construct local list */
    for ( int_t timestep=min_iteration; timestep<max_iteration ; timestep++ )
    {
        // Reinitialize list of actuals
        n_field = n_particles();
        resize_list ( n_field );

        // Prepare list of local particles
        marshal_particles ( &(list[0]) );

        // Add ghosts, append to list
        // This calculates n_virt, so n_field+n_virt=n_total for now
        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();
        generate_virtual_particles();
        t_end = MPI_Wtime();
        t_generate += t_end - t_start;

        // Synchronize with neighbors: mirror particles & mirror ghosts
        // Outcome is flat list, mirror particle count in n_mirror
        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();
        border_exchange();
        t_end = MPI_Wtime();
        t_border += t_end - t_start;

        // Node-local physics
        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();
        time_step ( timestep );
        t_end = MPI_Wtime();
        t_timestep += t_end - t_start;

        // Retrieve updates to actual particles from local list into hash tab
        unmarshal_particles( list, n_field );

        // Migrate particles moved across subdomain boundaries
        MPI_Barrier(MPI_COMM_WORLD);
        t_start = MPI_Wtime();
        migrate_particles();
        t_end = MPI_Wtime();
        t_migrate += t_end - t_start;

        // Write field state to file every few iterations
#ifndef NO_IO
        if ((timestep % checkpoint_frequency) == 0) {
            MPI_Barrier(MPI_COMM_WORLD);
            t_start = MPI_Wtime();
            char filename[256];
            memset ( filename, 0, 256*sizeof(char) );
            sprintf ( filename, "plot/%.4ld.dat",
                timestep/checkpoint_frequency
            );
            if ( rank == 0 )
                printf ( "Output at step %ld, '%s'\n", timestep, filename );

            collect_checkpoint();
            write_checkpoint ( filename );
            t_end = MPI_Wtime();
            t_io += t_end - t_start;
        } else if (verbose && rank == 0)
            printf("Step %ld\n", timestep);
#endif //NO_IO

    }

    /* Print all timings */
    print_timing("Generate ghosts: %.4lf, ", "%.4lf, ", t_generate);
    print_timing("Border exchange: %.4lf, ", "%.4lf, ", t_border);
    print_timing("Time step: %.4lf, ", "%.4lf, ", t_timestep);
    print_timing("Particle migration: %.4lf, ", "%.4lf, ", t_migrate);
    print_timing("Input/output: %.4lf, ", "%.4lf, ", t_io);
    print_timing("Find neighbors: %.4lf, ", "%.4lf, ", t_find_neighbors);

    /* Print shared variables */
    if (rank == 0) {
        printf ( "%lld particles\n", n_global_field );
        printf ( "%d ranks, OMP_NUM_THREADS=", size );
        int sys_retval = system ( "echo ${OMP_NUM_THREADS}" );
        printf ( "\t(returned %d)\n", sys_retval );
    }
}


int
main ( int argc, char **argv )
{
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &size );
    options ( argc, argv );
    east = (rank + 1) % size;
    west = (rank + size - 1) % size;

    if ( !restart )
        initialize();
    else
        restart_checkpoint ( min_iteration );

    time_integration();
    finalize();

    MPI_Finalize();
}


void
generate_virtual_particles ( void )
{
    n_virt = 0;
    real_t boundary = 1.55*H;

    // No particle adds more than 5 ghosts, make sure we have space
    resize_list(n_field * 5);

    // Loop over all actual particles, detect boundary interactions
    #pragma omp parallel for shared(n_virt)
    for ( int_t k=0; k<n_field; k++ )
    {
        int_t current_n_virt;

        // Horizontal mirror left
        if ( X(k) < boundary )
        {
            // Make use of one more list slot as ghost-k
            #pragma omp atomic capture
            current_n_virt = n_virt++;

            int_t gk = n_field + current_n_virt;
            X(gk) = -X(k), VX(gk) = -VX(k);
            Y(gk) = Y(k),  VY(gk) = VY(k);
            P(gk) = P(k), RHO(gk) = RHO(k), M(gk) = M(k);
            TYPE(gk) = -2, HSML(gk) = H;
        }

        // Horizontal mirror right
        if ( X(k) > B-boundary )
        {
            #pragma omp atomic capture
            current_n_virt = n_virt++;

            int_t gk = n_field + current_n_virt;
            X(gk) = 2*B-X(k), VX(gk) = -VX(k);
            Y(gk) = Y(k),     VY(gk) = VY(k);
            P(gk) = P(k), RHO(gk) = RHO(k);
            M(gk) = M(k); // Neumann boundary
            TYPE(gk) = -2, HSML(gk) = H;
        }

        // Vertical mirror bottom
        if ( Y(k) < boundary )
        {
            #pragma omp atomic capture
            current_n_virt = n_virt++;

            int_t gk = n_field + current_n_virt;
            X(gk) = X(k),  VX(gk) = VX(k);
            Y(gk) = -Y(k), VY(gk) = -VY(k);
            P(gk) = P(k), RHO(gk) = RHO(k), M(gk) = M(k);
            TYPE(gk) = -2, HSML(gk) = H;
        }

        // Lower left corner
        if ( X(k) < boundary && Y(k) < boundary )
        {
            #pragma omp atomic capture
            current_n_virt = n_virt++;

            int_t gk = n_field + current_n_virt;
            X(gk) = -X(k), VX(gk) = -VX(k);
            Y(gk) = -Y(k), VY(gk) = -VY(k);
            P(gk) = P(k), RHO(gk) = RHO(k), M(gk) = M(k);
            TYPE(gk) = -2, HSML(gk) = H;
        }

        // Lower right corner
        if ( X(k) > B-boundary && Y(k) < boundary )
        {
            #pragma omp atomic capture
            current_n_virt = n_virt++;

            int_t gk = n_field + current_n_virt;
            X(gk) = 2*B-X(k), VX(gk) = -VX(k);
            Y(gk) = -Y(k),    VY(gk) = -VY(k);
            P(gk) = P(k), RHO(gk) = RHO(k), M(gk) = M(k);
            TYPE(gk) = -2, HSML(gk) = H;
        }
    }
}


void
initialize ( void )
{
    /* Set up the hash table of actual (field) particles */
    particles_init ();

    /* Calculate the bounds of the local subdomain */
    real_t subdomain_size = B / (real_t)size;
    subdomain[0] = rank * subdomain_size;
    subdomain[1] = (rank+1) * subdomain_size;

    if ( subdomain_size <= sqrt((scale_k*H)*(scale_k*H)) )
        fprintf ( stderr, "Rank %d: "
            "Warning, subdomain size smaller than interaction radius at "
            "scale %lf, results will not be correct.\n", rank, SCALE
        );

    /* Populate the local subdomain with any initial particles */
    int_t n[2] = { 1+L/DELTA, 1+T/DELTA };
    n_global_field = n[0]*n[1];

    for ( int_t i=0; i<n[1]; i++ )
    {
        for ( int_t j=0; j<n[0]; j++ )
        {
            int_t k = j * n[1] + i; // k is particle index
            real_t
                x = H + j * DELTA,
                y = H + i * DELTA;
            if ( x >= subdomain[0] && x < subdomain[1] )
            {
                particle_t *p = malloc(sizeof(particle_t));
                p->idx = k;
                p->x[0] = x;
                p->x[1] = y;
                p->v[0] = p->v[1] = 0.0;
                p->rho = density;
                p->p = density * 9.81 * (T - p->x[1]);  // Hydrostatic pressure
                p->mass = L * T * density / (real_t)(n[0]*n[1]);
                p->type = 2;
                p->hsml = H;
                insert_particle ( p );
            }
        }
    }

    // Count local particles after initializeation
    n_field = n_particles();

    // Initial allocation for the lists of local particles and pairs
    list = (particle_t *) malloc ( n_capacity * sizeof(particle_t) );
    // Create one pair for every bucket
    pairs = (pair_t*)malloc(sizeof(pair_t) * N_BUCKETS_X * N_BUCKETS_Y);

    /* Create buckets */
    buckets = (bucket_t**)malloc(sizeof(bucket_t*) * N_BUCKETS_X * N_BUCKETS_Y);
    for (int i = 0; i < N_BUCKETS_X*N_BUCKETS_Y; ++i) {
        buckets[i] = (bucket_t*)malloc(sizeof(bucket_t));
    }


}


void
finalize ( void )
{
    particles_finalize ();
    free ( list );
    free ( pairs );
    free(buckets);

}


/* MPI communication */


void
migrate_particles ( void )
{
    particle_t
        *actuals[n_field],
        *move_west[n_field], // Can't migrate more than whole set
        *move_east[n_field];
    int_t export_east = 0, export_west = 0, import_east = 0, import_west = 0;

    list_particles ( actuals );
    for ( int_t k=0; k<n_field; k++ )
    {
        particle_t *p = actuals[k];
        if ( p->x[0] < subdomain[0] && rank > 0 )
        {
            remove_particle ( p );
            move_west[export_west] = p;
            export_west += 1;
        }
        else if ( p->x[0] > subdomain[1] && rank < (size-1) )
        {
            remove_particle ( p );
            move_east[export_east] = p;
            export_east += 1;
        }
/*
        if ( p->x[0] < 0.0 )
            p->x[0] += B;
        else if ( p->x[0] > B )
            p->x[0] -= B;
*/
    }
    MPI_Sendrecv (
        &export_west, 1, INT_MACRO_MPI, west, 0,
        &import_east, 1, INT_MACRO_MPI, east, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv (
        &export_east, 1, INT_MACRO_MPI, east, 0,
        &import_west, 1, INT_MACRO_MPI, west, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );

    // Serialize in/out-bound particles
    particle_t *outlist = malloc (
        (export_west+export_east)*sizeof(particle_t)
    );
    particle_t *inlist = malloc (
        (import_west+import_east)*sizeof(particle_t)
    );

    // Outbound
    for ( int_t k=0; k<export_west; k++ )
        memcpy ( &(outlist[k]), move_west[k], sizeof(particle_t) );
    for ( int_t k=0; k<export_east; k++ )
        memcpy ( &(outlist[export_west+k]), move_east[k], sizeof(particle_t) );

    MPI_Sendrecv (
        &(outlist[0]),
        export_west*sizeof(particle_t), MPI_BYTE, west, 0,
        &(inlist[import_west]),
        import_east*sizeof(particle_t), MPI_BYTE, east, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv (
        &(outlist[export_west]),
        export_east*sizeof(particle_t), MPI_BYTE, east, 0,
        &(inlist[0]),
        import_west*sizeof(particle_t), MPI_BYTE, west, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );

    for ( int_t k=0; k<(import_west+import_east); k++ )
    {
        particle_t *p = malloc ( sizeof(particle_t) );
        memcpy ( p, &(inlist[k]), sizeof(particle_t) );
        insert_particle ( p );
    }

    free ( outlist );
    free ( inlist );

    // Remove migrated particles from local responsibility
    for ( int_t k=0; k<export_west; k++ )
        free ( move_west[k] );
    for ( int_t k=0; k<export_east; k++ )
        free ( move_east[k] );
}


void
border_exchange ( void )
{
    int_t
        export_east = 0, export_west = 0,
        import_east = 0, import_west = 0;

    // Scan the list and count particles within neighbor reach
    #pragma omp parallel for
    for ( int_t k=0; k<(n_field+n_virt); k++ )
    {
        if ( (X(k) - subdomain[0]) < scale_k*H && rank>0 )
            #pragma omp atomic update
            export_west += 1;
        if ( (subdomain[1] - X(k)) < scale_k*H && rank<size-1)
            #pragma omp atomic update
            export_east += 1;
    }

    MPI_Sendrecv (
        &export_west, 1, INT_MACRO_MPI, west, 0,
        &import_east, 1, INT_MACRO_MPI, east, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv (
        &export_east, 1, INT_MACRO_MPI, east, 0,
        &import_west, 1, INT_MACRO_MPI, west, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );

    // This transfer list could be glob/resize instead of malloc per iter

    particle_t *transfer = (particle_t *) malloc (
        (export_west + export_east) * sizeof(particle_t)
    );
    int_t w_idx = 0, e_idx=export_west;
    int_t priv_w_idx,  priv_e_idx;

    #pragma omp parallel for private(priv_w_idx, priv_e_idx)
    for ( int_t k=0; k<(n_field+n_virt); k++ )
    {
        if ( (X(k) - subdomain[0]) < scale_k*H && rank > 0)
        {
            #pragma omp atomic capture
            priv_w_idx = w_idx++;
            memcpy ( &(transfer[priv_w_idx]), &(list[k]), sizeof(particle_t) );
        }
        if ( (subdomain[1] - X(k)) < scale_k*H && rank < size-1)
        {
            #pragma omp atomic capture
            priv_e_idx = e_idx++;
            memcpy ( &(transfer[priv_e_idx]), &(list[k]), sizeof(particle_t) );
        }
    }

    n_mirror = import_east + import_west;
    resize_list ( n_field + n_virt + n_mirror );

    MPI_Sendrecv (
        &(transfer[0]),
        export_west*sizeof(particle_t), MPI_BYTE, west, 0,
        &(list[n_field+n_virt+import_west]),
        import_east*sizeof(particle_t), MPI_BYTE, east, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv (
        &(transfer[export_west]),
        export_east*sizeof(particle_t), MPI_BYTE, east, 0,
        &(list[n_field+n_virt]),
        import_west*sizeof(particle_t), MPI_BYTE, west, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE
    );
    free ( transfer );
}

/* Auxiliary routines - file handling is in sph_io.c */

void
resize_list ( int_t required )
{
    if ( required >= n_capacity )
    {
        n_capacity = required;
        particle_t* new_list = realloc ( list, n_capacity * sizeof(particle_t) );
        if(new_list == NULL) {
            free(list);
            fprintf(stderr, "1Not enough memory!\n");
            exit(1);
        } else {
            list = new_list;
        }
    }
}


void
resize_pair_list ( int_t new_cap )
{
    n_pair_cap = new_cap;
    pair_t* new_pairs = realloc ( pairs, n_pair_cap * sizeof(pair_t) );
    if(new_pairs == NULL) {
        fprintf(stderr, "2Not enough memory!\n");
        free(pairs);
        exit(1);
    } else {
        pairs = new_pairs;
    }
}


void
options ( int argc, char **argv )
{
    /* Master rank parses command line options */
    if ( rank == 0 )
    {
        int o;
        while ( (o = getopt(argc,argv,"i:c:r:")) != -1 )
        switch ( o )
        {
            case 'i':
                max_iteration =
                    strtol(optarg, NULL, 10);
                break;
            case 'c':
                checkpoint_frequency =
                    strtol(optarg, NULL, 10);
                break;
            case 'r':
                min_iteration = strtol(optarg,NULL,10);
                break;
        }
    }

    /* Communicate option flags to the rest of the collective */
    MPI_Bcast ( &min_iteration, 1, INT_MACRO_MPI, 0,
        MPI_COMM_WORLD
    );
    MPI_Bcast ( &max_iteration, 1, INT_MACRO_MPI, 0,
        MPI_COMM_WORLD
    );
    MPI_Bcast ( &checkpoint_frequency, 1, INT_MACRO_MPI, 0,
        MPI_COMM_WORLD
    );
    if ( min_iteration != MIN_ITERATION_DEFAULT )
        restart = true;
}
