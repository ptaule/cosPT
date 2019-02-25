using Test, DelimitedFiles

const vfloat = Float64

const datasets_dir = "/home/t30/ben/ge52sir/Code/Non_linear_PS/tests/mathematica/"

const points = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10];

const mu_vals = [0,0.3,0.5,0.8,0.9];
filenames = ["F2_mu_0.dat", "F2_mu_03.dat", "F2_mu_05.dat", "F2_mu_08.dat", "F2_mu_09.dat" ];

const abs_tol = 1e-10
const rel_tol = 1e-10

@testset "F2kernel" begin
    n = 2;
    component = 0;
    filenames = ["F2_mu_0.dat", "F2_mu_03.dat", "F2_mu_05.dat", "F2_mu_08.dat", "F2_mu_09.dat" ];

    for a=1:5
        μ = mu_vals[a]
        mathematica = readdlm(datasets_dir * "kernels/" * filenames[a]);

        @testset "μ = $μ" begin
            for b=1:6, c=1:6
                k = points[b];
                q = points[c];
                # Call c-interface to kernel-computer
                kernel = ccall((:SPTkernel1Loop,
                                "/home/t30/ben/ge52sir/Code/Non_linear_PS/tests/test_interface.so"),
                               vfloat,(Int32,Int32,vfloat,vfloat,vfloat),n,component,k,q,μ);
                @test isapprox(kernel,mathematica[b,c];atol=abs_tol,rtol=rel_tol)
            end
        end
    end
end

@testset "G2kernel" begin
    n = 2;
    component = 1;

    filenames = ["G2_mu_0.dat", "G2_mu_03.dat", "G2_mu_05.dat", "G2_mu_08.dat", "G2_mu_09.dat" ];
    for a=1:5
        μ = mu_vals[a]
        mathematica = readdlm(datasets_dir * "kernels/" * filenames[a]);

        @testset "μ = $μ" begin
            for b=1:6, c=1:6
                k = points[b];
                q = points[c];
                # Call c-interface to kernel-computer
                kernel = ccall((:SPTkernel1Loop,
                                "/home/t30/ben/ge52sir/Code/Non_linear_PS/tests/test_interface.so"),
                               vfloat,(Int32,Int32,vfloat,vfloat,vfloat),n,component,k,q,μ);

                @test isapprox(kernel,mathematica[b,c];atol=abs_tol,rtol=rel_tol)
            end
        end
    end
end

@testset "F3kernel" begin
    n = 3;
    component = 0;

    filenames = ["F3_mu_0.dat", "F3_mu_03.dat", "F3_mu_05.dat", "F3_mu_08.dat", "F3_mu_09.dat" ];
    for a=1:5
        μ = mu_vals[a]
        mathematica = readdlm(datasets_dir * "kernels/" * filenames[a]);

        @testset "μ = $μ" begin
            for b=1:6, c=1:6
                k = points[b];
                q = points[c];
                # Call c-interface to kernel-computer
                kernel = ccall((:SPTkernel1Loop,
                                "/home/t30/ben/ge52sir/Code/Non_linear_PS/tests/test_interface.so"),
                               vfloat,(Int32,Int32,vfloat,vfloat,vfloat),n,component,k,q,μ);

                @test isapprox(kernel,mathematica[b,c];atol=abs_tol,rtol=rel_tol)
            end
        end
    end
end

@testset "G3kernel" begin
    n = 3;
    component = 1;

    filenames = ["G3_mu_0.dat", "G3_mu_03.dat", "G3_mu_05.dat", "G3_mu_08.dat", "G3_mu_09.dat" ];
    for a=1:5
        μ = mu_vals[a]
        mathematica = readdlm(datasets_dir * "kernels/" * filenames[a]);

        @testset "μ = $μ" begin
            for b=1:6, c=1:6
                k = points[b];
                q = points[c];
                # Call c-interface to kernel-computer
                kernel = ccall((:SPTkernel1Loop,
                                "/home/t30/ben/ge52sir/Code/Non_linear_PS/tests/test_interface.so"),
                               vfloat,(Int32,Int32,vfloat,vfloat,vfloat),n,component,k,q,μ);

                @test isapprox(kernel,mathematica[b,c];atol=abs_tol,rtol=rel_tol)
            end
        end
    end
end
