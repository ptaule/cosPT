


void testIntegrandComputation() {
    gsl_interp_accel* acc;
    gsl_spline* spline;
    const char* filename = "/home/pettertaule/Dropbox/Mathematica/simple00_pk.dat";

    read_PS(filename,&acc,&spline);

    integration_variables_t vars = {
        .magnitudes = {0.2},
        .cos_theta  = {1},
        .phi        = {}
    };

    integration_input_t data = {
        .k = 1,
        .component_a = 0,
        .component_b = 0,
        .acc = acc,
        .spline = spline
    };

    vfloat result = integrand(&data,&vars);
    printf("result  = %f\n", result );

}



void testKernelComputer() {
    vfloat k = 3;
    integration_variables_t vars = {
        .magnitudes = {1,2},
        .cos_theta  = {0.1,0.5},
        .phi        = {0}
    };

    short int component = 0;
    short int n = 5;

    short int args[N_KERNEL_ARGS] = {13,5,3,7,1};
    // 2-loop:
    // k       <-> 13
    // k - Q_2 <-> 10
    // Q_1     <-> 5
    // - Q_1   <-> 3
    // Q_2     <-> 7
    // - Q_2   <-> 1

    table_pointers_t data_tables;
    data_tables.Q_magnitudes = vars.magnitudes,
    data_tables.alpha = matrix_alloc(N_CONFIGS,N_CONFIGS);
    data_tables.beta  = matrix_alloc(N_CONFIGS,N_CONFIGS);

    // Allocate space for kernels (calloc also initializes values to 0)
    data_tables.kernels = (kernel_value_t*)
        calloc(COMPONENTS * N_KERNELS, sizeof(kernel_value_t));

    // Initialize sum-, bare_scalar_products-, alpha- and beta-tables
    compute_sum_table(data_tables.sum_table);
    compute_bare_scalar_products(k, &vars,
            data_tables.bare_scalar_products);
    compute_alpha_beta_tables(data_tables.bare_scalar_products,
            data_tables.alpha, data_tables.beta);

    vfloat value = compute_SPT_kernel(args,n,component,&data_tables);
    printf("result = " vfloat_fmt "\n",value);

    // Free allocated memory
    free(data_tables.kernels);
    matrix_free(data_tables.alpha);
    matrix_free(data_tables.beta);
}



void print_configs() {
    short int configs[N_COEFFS];

    for (short int label = 0; label < N_CONFIGS; ++label) {
        label2config(label,configs,N_COEFFS);
        for (int i = 0; i < N_COEFFS; ++i) {
            printf("%d,",configs[i]);
        }
        printf("\n");
    }
}


void testAlphaBeta() {
    vfloat k = 2;
    integration_variables_t vars = {
        .magnitudes = {1,2},
        .cos_theta  = {0.1,0.5},
        .phi        = {0}
    };

    vfloat bare_scalar_products[N_COEFFS][N_COEFFS];
    matrix_t* alpha = matrix_alloc(N_CONFIGS,N_CONFIGS);
    matrix_t* beta =  matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_bare_scalar_products(k,&vars,bare_scalar_products);
    compute_alpha_beta_tables(bare_scalar_products,alpha,beta);

    printf("alpha=\n");
    print_gsl_matrix(alpha,N_CONFIGS,N_CONFIGS);
    printf("beta=\n");
    print_gsl_matrix(beta,N_CONFIGS,N_CONFIGS);

    matrix_free(alpha);
    matrix_free(beta);
}



void testVectorSum() {
    short int args[3] = {9,7,5};

    for (int j = 0; j < 3; ++j) {
        short int config[N_COEFFS];
        label2config(args[j],config,N_COEFFS);
        for (int i = 0; i < N_COEFFS; ++i) {
            printf("%d,",config[i]);
        }
        printf(" + ");
    }
    printf("=\n");

    short int sum = sum_vectors(args,3);

    short int sumConfig[N_COEFFS];
    label2config(sum,sumConfig,N_COEFFS);

    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",sumConfig[i]);
    }
    printf("\n");
}




void test_kernel_index_from_arguments() {
    short int label1 = 17;
    short int label2 = 5;
    short int label3 = 1;
    short int label4 = 3;
    short int label5 = 7;

    short int arguments[5] = {label1,label2,label3,label4,label5};
    short int config1[N_COEFFS];
    short int config2[N_COEFFS];
    short int config3[N_COEFFS];
    short int config4[N_COEFFS];
    short int config5[N_COEFFS];

    label2config(label1,config1,N_COEFFS);
    label2config(label2,config2,N_COEFFS);
    label2config(label3,config3,N_COEFFS);
    label2config(label4,config4,N_COEFFS);
    label2config(label5,config5,N_COEFFS);

    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config1[i]);
    }
    printf("\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config2[i]);
    }
    printf("\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config3[i]);
    }
    printf("\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config4[i]);
    }
    printf("\n");
    for (int i = 0; i < N_COEFFS; ++i) {
        printf("%d,",config5[i]);
    }
    printf("\n");

    short int index = 0;
    short int n = 0;
    kernel_index_from_arguments(arguments,&index,&n);

    printf("kernel_index_from_arguments = %d\n",index);
}



void testBareScalarProducts() {
    vfloat bare_scalar_products[N_COEFFS][N_COEFFS] = {};

    vfloat k = 2;
    integration_variables_t vars = {
        .magnitudes = {1,2},
        .cos_theta  = {0.1,0.5},
        .phi        = {0}
    };

    compute_bare_scalar_products(k,&vars,bare_scalar_products);

    for (int i = 0; i < N_COEFFS; ++i) {
        for (int j = 0; j < N_COEFFS; ++j) {
            printf("%2.5f  ",bare_scalar_products[i][j]);
        }
        printf("\n");
    }
}



void testScalarProducts() {
    vfloat k = 2;
    integration_variables_t vars = {
        .magnitudes = {1,2},
        .cos_theta  = {0.1,0.5},
        .phi        = {0}
    };

    matrix_t* scalar_products = matrix_alloc(N_CONFIGS,N_CONFIGS);
    vfloat bare_scalar_products[N_COEFFS][N_COEFFS];

    compute_bare_scalar_products(k,&vars,bare_scalar_products);
    compute_scalar_products(bare_scalar_products,scalar_products);

    for (int i = 0; i < N_CONFIGS; ++i) {
        for (int j = 0; j < N_CONFIGS; ++j) {
            printf("%2.2f  ",matrix_get(scalar_products,i,j));
        }
        printf("\n");
    }
}
